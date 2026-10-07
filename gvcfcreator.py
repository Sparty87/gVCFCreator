#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# gVCFCreator
# Copyright (c) 2026 Serena Spartano. All rights reserved.
"""
gVCFCreator
Generation of genomic VCFs (g.vcf.gz) from BAM/BAI files using GATK HaplotypeCaller (-ERC GVCF).

Requirements:
  - Python 3.8+ with Tkinter
  - GATK 4.x ('gatk' launcher or .jar file) + Java 17
  - samtools (optional, only needed to create the reference .fai if missing)
"""

import os
import re
import gzip
import struct
import shutil
import shlex
import subprocess
import threading
import queue
import time
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter.scrolledtext import ScrolledText

__author__ = "Serena Spartano"
__copyright__ = "Copyright (c) 2026 Serena Spartano"
__version__ = "1.0.0"

APP_TITLE = "gVCFCreator — BAM → gVCF (GATK HaplotypeCaller)"


# ----------------------------------------------------------------------------- utilities
def find_bai(bam):
    """Return the path of the BAM index if it exists (sample.bam.bai or sample.bai)."""
    for cand in (bam + ".bai", os.path.splitext(bam)[0] + ".bai", bam[:-4] + ".csi" if bam.endswith(".bam") else None):
        if cand and os.path.isfile(cand):
            return cand
    return None


def ref_dict_path(ref):
    base = ref
    for ext in (".gz", ".fasta", ".fa", ".fna"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base + ".dict"


def sample_name(bam):
    name = os.path.basename(bam)
    exts = (".bam", ".cram", ".sorted", ".dedup", ".markdup", ".recal", ".bqsr")
    changed = True
    while changed:
        changed = False
        for ext in exts:
            if name.endswith(ext) and len(name) > len(ext):
                name = name[: -len(ext)]
                changed = True
    return name


def read_bam_header(path):
    """Read the text header of a BAM (BGZF) with no external dependencies; for CRAM, use samtools."""
    if path.endswith(".bam"):
        with gzip.open(path, "rb") as f:
            if f.read(4) != b"BAM\x01":
                raise ValueError("File not recognized as BAM")
            l_text = struct.unpack("<i", f.read(4))[0]
            return f.read(l_text).decode("utf-8", errors="replace")
    if not shutil.which("samtools"):
        raise RuntimeError("samtools is required to read the CRAM header")
    return subprocess.run(["samtools", "view", "-H", path], capture_output=True,
                          text=True, check=True).stdout


def bam_samples(path):
    """Return the list of samples (SM tag of the @RG read groups), in order of appearance."""
    samples = []
    for line in read_bam_header(path).splitlines():
        if line.startswith("@RG"):
            for field in line.split("\t")[1:]:
                if field.startswith("SM:") and field[3:] not in samples:
                    samples.append(field[3:])
    return samples


def safe_name(txt):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", txt).strip("_") or "sample"


# ----------------------------------------------------------------------------- GUI
class GVCFApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("980x760")
        self.minsize(860, 640)

        self.msg_q = queue.Queue()
        self.proc = None
        self.worker = None
        self.stop_flag = threading.Event()

        # variables
        self.v_gatk = tk.StringVar(value=shutil.which("gatk") or "")
        self.v_mem = tk.IntVar(value=8)
        self.v_ref = tk.StringVar()
        self.v_out = tk.StringVar()
        self.v_intervals = tk.StringVar()
        self.v_padding = tk.IntVar(value=0)
        self.v_mapq = tk.IntVar(value=20)
        self.v_bq = tk.IntVar(value=10)
        self.v_conf = tk.DoubleVar(value=30.0)
        self.v_ploidy = tk.IntVar(value=2)
        self.v_threads = tk.IntVar(value=4)
        self.v_dups = tk.BooleanVar(value=True)
        self.v_make_idx = tk.BooleanVar(value=True)
        self.v_overwrite = tk.BooleanVar(value=False)
        self.v_multi = tk.StringVar(value="split")

        self._build()
        self.after(150, self._poll_queue)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    # ------------------------------------------------------------------ layout
    def _build(self):
        pad = {"padx": 6, "pady": 4}
        main = ttk.Frame(self, padding=8)
        main.pack(fill="both", expand=True)

        # --- paths
        f_paths = ttk.LabelFrame(main, text="Paths")
        f_paths.pack(fill="x", **pad)
        f_paths.columnconfigure(1, weight=1)

        self._path_row(f_paths, 0, "GATK (launcher or .jar):", self.v_gatk, self._browse_gatk)
        ttk.Label(f_paths, text="Java memory (GB):").grid(row=0, column=3, sticky="e", **pad)
        ttk.Spinbox(f_paths, from_=2, to=256, textvariable=self.v_mem, width=6).grid(row=0, column=4, **pad)

        self._path_row(f_paths, 1, "Reference genome (FASTA):", self.v_ref, self._browse_ref)
        ttk.Checkbutton(f_paths, text="Create .fai/.dict if missing",
                        variable=self.v_make_idx).grid(row=1, column=3, columnspan=2, sticky="w", **pad)

        self._path_row(f_paths, 2, "Output folder:", self.v_out, self._browse_out)
        ttk.Checkbutton(f_paths, text="Overwrite existing files",
                        variable=self.v_overwrite).grid(row=2, column=3, columnspan=2, sticky="w", **pad)

        self._path_row(f_paths, 3, "Regions (BED/interval_list, optional):", self.v_intervals, self._browse_intervals)
        ttk.Label(f_paths, text="Padding (bp):").grid(row=3, column=3, sticky="e", **pad)
        ttk.Spinbox(f_paths, from_=0, to=1000, textvariable=self.v_padding, width=6).grid(row=3, column=4, **pad)

        # --- BAM
        f_bam = ttk.LabelFrame(main, text="BAM files (the .bai index must be located next to the BAM)")
        f_bam.pack(fill="x", **pad)
        f_bam.columnconfigure(0, weight=1)
        self.lb_bam = tk.Listbox(f_bam, height=6, selectmode="extended")
        self.lb_bam.grid(row=0, column=0, rowspan=3, sticky="nsew", **pad)
        sb = ttk.Scrollbar(f_bam, orient="vertical", command=self.lb_bam.yview)
        sb.grid(row=0, column=1, rowspan=3, sticky="ns")
        self.lb_bam.config(yscrollcommand=sb.set)
        ttk.Button(f_bam, text="Add BAM…", command=self._add_bams).grid(row=0, column=2, sticky="ew", **pad)
        ttk.Button(f_bam, text="Remove selected", command=self._remove_bams).grid(row=1, column=2, sticky="ew", **pad)
        ttk.Button(f_bam, text="Clear", command=lambda: self.lb_bam.delete(0, "end")).grid(row=2, column=2, sticky="ew", **pad)

        # --- quality filters
        f_q = ttk.LabelFrame(main, text="Quality filters and parameters")
        f_q.pack(fill="x", **pad)
        params = [
            ("Minimum read MAPQ:", self.v_mapq, 0, 255, 1,
             "--minimum-mapping-quality: discards reads with lower MAPQ (GATK default 20)"),
            ("Minimum base quality:", self.v_bq, 0, 60, 1,
             "--min-base-quality-score: bases below threshold are ignored (GATK default 10)"),
            ("Calling confidence threshold (Phred):", self.v_conf, 0, 100, 1,
             "--standard-min-confidence-threshold-for-calling (default 30)"),
            ("Ploidy:", self.v_ploidy, 1, 10, 1, "--sample-ploidy (2 for human autosomes)"),
            ("PairHMM threads:", self.v_threads, 1, 64, 1, "--native-pair-hmm-threads"),
        ]
        for i, (lab, var, lo, hi, inc, tip) in enumerate(params):
            r, c = divmod(i, 2)
            ttk.Label(f_q, text=lab).grid(row=r, column=c * 3, sticky="e", **pad)
            ttk.Spinbox(f_q, from_=lo, to=hi, increment=inc, textvariable=var, width=7).grid(row=r, column=c * 3 + 1, sticky="w", **pad)
            ttk.Label(f_q, text=tip, foreground="gray40", font=("TkDefaultFont", 8)).grid(row=r, column=c * 3 + 2, sticky="w", **pad)
        ttk.Checkbutton(f_q, text="Exclude reads marked as duplicates (GATK default filter)",
                        variable=self.v_dups).grid(row=3, column=0, columnspan=6, sticky="w", **pad)
        f_ms = ttk.Frame(f_q)
        f_ms.grid(row=4, column=0, columnspan=6, sticky="w", **pad)
        ttk.Label(f_ms, text="BAM with multiple samples (different SM tags):").pack(side="left")
        ttk.Radiobutton(f_ms, text="one gVCF per sample (--sample-name)",
                        variable=self.v_multi, value="split").pack(side="left", padx=6)
        ttk.Radiobutton(f_ms, text="merge into a single sample (AddOrReplaceReadGroups)",
                        variable=self.v_multi, value="merge").pack(side="left", padx=6)

        # --- buttons
        f_btn = ttk.Frame(main)
        f_btn.pack(fill="x", **pad)
        ttk.Button(f_btn, text="Validate inputs", command=self._validate_only).pack(side="left", padx=4)
        self.btn_run = ttk.Button(f_btn, text="▶ Run", command=self._start)
        self.btn_run.pack(side="left", padx=4)
        self.btn_stop = ttk.Button(f_btn, text="■ Stop", command=self._stop, state="disabled")
        self.btn_stop.pack(side="left", padx=4)
        self.pb = ttk.Progressbar(f_btn, mode="determinate")
        self.pb.pack(side="left", fill="x", expand=True, padx=8)
        self.lbl_status = ttk.Label(f_btn, text="Ready")
        self.lbl_status.pack(side="right", padx=4)

        # --- log
        f_log = ttk.LabelFrame(main, text="Log")
        f_log.pack(fill="both", expand=True, **pad)
        self.txt = ScrolledText(f_log, height=14, font=("Courier", 9), wrap="none")
        self.txt.pack(fill="both", expand=True)

    def _path_row(self, parent, row, label, var, cmd):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="e", padx=6, pady=4)
        ttk.Entry(parent, textvariable=var).grid(row=row, column=1, sticky="ew", padx=6, pady=4)
        ttk.Button(parent, text="Browse…", command=cmd).grid(row=row, column=2, padx=6, pady=4)

    # ------------------------------------------------------------------ browse
    def _browse_gatk(self):
        p = filedialog.askopenfilename(title="Select gatk or gatk-package-*.jar")
        if p:
            self.v_gatk.set(p)

    def _browse_ref(self):
        p = filedialog.askopenfilename(title="Select reference FASTA",
                                       filetypes=[("FASTA", "*.fa *.fasta *.fna *.fa.gz *.fasta.gz"), ("All files", "*")])
        if p:
            self.v_ref.set(p)

    def _browse_out(self):
        p = filedialog.askdirectory(title="Select output folder")
        if p:
            self.v_out.set(p)

    def _browse_intervals(self):
        p = filedialog.askopenfilename(title="Select regions file",
                                       filetypes=[("Intervals", "*.bed *.interval_list *.intervals *.list"), ("All files", "*")])
        if p:
            self.v_intervals.set(p)

    def _add_bams(self):
        files = filedialog.askopenfilenames(title="Select BAM", filetypes=[("BAM", "*.bam"), ("CRAM", "*.cram")])
        current = set(self.lb_bam.get(0, "end"))
        for f in files:
            if f not in current:
                self.lb_bam.insert("end", f)

    def _remove_bams(self):
        for i in reversed(self.lb_bam.curselection()):
            self.lb_bam.delete(i)

    # ------------------------------------------------------------------ log
    def log(self, msg):
        self.msg_q.put(("log", msg))

    def _poll_queue(self):
        try:
            while True:
                kind, payload = self.msg_q.get_nowait()
                if kind == "log":
                    self.txt.insert("end", payload.rstrip("\n") + "\n")
                    self.txt.see("end")
                elif kind == "progress":
                    done, tot, label = payload
                    self.pb.config(maximum=max(tot, 1), value=done)
                    self.lbl_status.config(text=label)
                elif kind == "done":
                    self.btn_run.config(state="normal")
                    self.btn_stop.config(state="disabled")
                    self.lbl_status.config(text=payload)
        except queue.Empty:
            pass
        self.after(150, self._poll_queue)

    # ------------------------------------------------------------------ validation
    def _gatk_base_cmd(self):
        gatk = self.v_gatk.get().strip()
        mem = self.v_mem.get()
        if gatk.endswith(".jar"):
            return ["java", f"-Xmx{mem}g", "-jar", gatk]
        return [gatk, "--java-options", f"-Xmx{mem}g"]

    def _collect_errors(self):
        errs, warns = [], []
        gatk = self.v_gatk.get().strip()
        if not gatk:
            errs.append("GATK path not specified.")
        elif not (os.path.isfile(gatk) or shutil.which(gatk)):
            errs.append(f"GATK not found: {gatk}")
        elif gatk.endswith(".jar") and not shutil.which("java"):
            errs.append("Java not found in PATH (required to run the .jar).")

        ref = self.v_ref.get().strip()
        if not ref or not os.path.isfile(ref):
            errs.append("Invalid reference genome.")
        else:
            if ref.endswith(".gz") and not os.path.isfile(ref + ".gzi"):
                warns.append("Compressed FASTA without .gzi: GATK requires bgzip + .gzi index.")
            if not os.path.isfile(ref + ".fai"):
                (warns if self.v_make_idx.get() else errs).append(f"Missing .fai index: {ref}.fai")
            if not os.path.isfile(ref_dict_path(ref)):
                (warns if self.v_make_idx.get() else errs).append(f"Missing sequence dictionary: {ref_dict_path(ref)}")
            if self.v_make_idx.get() and not os.path.isfile(ref + ".fai") and not shutil.which("samtools"):
                errs.append("samtools not found: unable to create the .fai.")

        out = self.v_out.get().strip()
        if not out:
            errs.append("Output folder not specified.")

        itv = self.v_intervals.get().strip()
        if itv and not os.path.isfile(itv):
            errs.append(f"Regions file not found: {itv}")

        bams = list(self.lb_bam.get(0, "end"))
        if not bams:
            errs.append("No BAM file selected.")
        for b in bams:
            if not os.path.isfile(b):
                errs.append(f"BAM not found: {b}")
            elif b.endswith(".bam") and not find_bai(b):
                errs.append(f"Missing .bai index for: {os.path.basename(b)}")
            elif b.endswith(".cram") and not os.path.isfile(b + ".crai"):
                errs.append(f"Missing .crai index for: {os.path.basename(b)}")
        names = [sample_name(b) for b in bams]
        if len(set(names)) != len(names):
            errs.append("Multiple BAMs produce the same output name: please rename the files.")
        return errs, warns

    def _validate_only(self):
        errs, warns = self._collect_errors()
        for w in warns:
            self.log(f"[WARNING] {w}")
        for e in errs:
            self.log(f"[ERROR] {e}")
        if errs:
            messagebox.showerror("Validation", f"{len(errs)} error(s). See the log.")
        else:
            try:
                out = subprocess.run(self._gatk_base_cmd() + ["--version"], capture_output=True, text=True, timeout=120)
                self.log((out.stdout + out.stderr).strip())
            except Exception as ex:
                self.log(f"[WARNING] Unable to query the GATK version: {ex}")
            messagebox.showinfo("Validation", "Inputs are valid." + (f"\n{len(warns)} warning(s) in the log." if warns else ""))

    # ------------------------------------------------------------------ execution
    def _start(self):
        errs, warns = self._collect_errors()
        for w in warns:
            self.log(f"[WARNING] {w}")
        if errs:
            for e in errs:
                self.log(f"[ERROR] {e}")
            messagebox.showerror("Error", "Please fix the errors reported in the log.")
            return
        os.makedirs(self.v_out.get().strip(), exist_ok=True)
        cfg = {
            "base": self._gatk_base_cmd(),
            "ref": self.v_ref.get().strip(),
            "out": self.v_out.get().strip(),
            "itv": self.v_intervals.get().strip(),
            "pad": self.v_padding.get(),
            "mapq": self.v_mapq.get(),
            "bq": self.v_bq.get(),
            "conf": self.v_conf.get(),
            "ploidy": self.v_ploidy.get(),
            "threads": self.v_threads.get(),
            "dups": self.v_dups.get(),
            "make_idx": self.v_make_idx.get(),
            "overwrite": self.v_overwrite.get(),
            "multi": self.v_multi.get(),
            "bams": list(self.lb_bam.get(0, "end")),
        }
        self.stop_flag.clear()
        self.btn_run.config(state="disabled")
        self.btn_stop.config(state="normal")
        self.worker = threading.Thread(target=self._run_all, args=(cfg,), daemon=True)
        self.worker.start()

    def _stop(self):
        self.stop_flag.set()
        if self.proc and self.proc.poll() is None:
            self.log("[INFO] Stop requested, terminating the running process…")
            self.proc.terminate()

    def _exec(self, cmd, logfile=None):
        self.log("$ " + " ".join(shlex.quote(c) for c in cmd))
        fh = open(logfile, "w") if logfile else None
        try:
            self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                         text=True, bufsize=1)
            for line in self.proc.stdout:
                if fh:
                    fh.write(line)
                # reduce noise: show only the relevant informative lines
                if any(k in line for k in ("ProgressMeter", "ERROR", "WARN", "Exception",
                                            "HaplotypeCaller done", "Traversal complete")):
                    self.log(line)
            return self.proc.wait()
        finally:
            if fh:
                fh.close()
            self.proc = None

    def _prepare_reference(self, cfg):
        ref = cfg["ref"]
        if not cfg["make_idx"]:
            return True
        if not os.path.isfile(ref + ".fai"):
            self.log("[INFO] Creating .fai index…")
            if self._exec(["samtools", "faidx", ref]) != 0:
                self.log("[ERROR] samtools faidx failed.")
                return False
        if not os.path.isfile(ref_dict_path(ref)):
            self.log("[INFO] Creating .dict sequence dictionary…")
            if self._exec(cfg["base"] + ["CreateSequenceDictionary", "-R", ref, "-O", ref_dict_path(ref)]) != 0:
                self.log("[ERROR] CreateSequenceDictionary failed.")
                return False
        return True

    def _build_hc_cmd(self, cfg, bam, out_gvcf, sm=None):
        cmd = cfg["base"] + [
            "HaplotypeCaller",
            "-R", cfg["ref"],
            "-I", bam,
            "-O", out_gvcf,
            "-ERC", "GVCF",
            "--minimum-mapping-quality", str(cfg["mapq"]),
            "--min-base-quality-score", str(cfg["bq"]),
            "--standard-min-confidence-threshold-for-calling", str(cfg["conf"]),
            "--sample-ploidy", str(cfg["ploidy"]),
            "--native-pair-hmm-threads", str(cfg["threads"]),
        ]
        if cfg["itv"]:
            cmd += ["-L", cfg["itv"]]
            if cfg["pad"] > 0:
                cmd += ["-ip", str(cfg["pad"])]
        if not cfg["dups"]:
            cmd += ["--disable-read-filter", "NotDuplicateReadFilter"]
        if sm:
            cmd += ["--sample-name", sm]
        return cmd

    def _unify_read_groups(self, cfg, bam, name):
        """Rewrite all read groups into a single one with SM=name (temporary BAM in the output folder)."""
        tmp_bam = os.path.join(cfg["out"], f".{name}.tmp_rg.bam")
        cmd = cfg["base"] + [
            "AddOrReplaceReadGroups", "-I", bam, "-O", tmp_bam,
            "--RGID", name, "--RGSM", name, "--RGLB", "lib1", "--RGPL", "ILLUMINA", "--RGPU", "unit1",
            "--CREATE_INDEX", "true", "--VALIDATION_STRINGENCY", "SILENT",
        ]
        self.log(f"[INFO] Merging read groups → SM={name}")
        rc = self._exec(cmd, os.path.join(cfg["out"], f"{name}.addrg.log"))
        return tmp_bam if rc == 0 and os.path.isfile(tmp_bam) else None

    @staticmethod
    def _remove_tmp(tmp_bam):
        if tmp_bam:
            for f in (tmp_bam, tmp_bam[:-4] + ".bai", tmp_bam + ".bai"):
                if os.path.isfile(f):
                    os.remove(f)

    def _run_all(self, cfg):
        t0 = time.time()
        bams = cfg["bams"]
        ok, fail, skip = [], [], []
        self.msg_q.put(("progress", (0, len(bams), "Preparing reference…")))

        if not self._prepare_reference(cfg):
            self.msg_q.put(("done", "Reference error"))
            return

        for i, bam in enumerate(bams):
            if self.stop_flag.is_set():
                break
            name = sample_name(bam)
            self.msg_q.put(("progress", (i, len(bams), f"[{i + 1}/{len(bams)}] {name}")))
            self.log(f"\n===== {name} ({i + 1}/{len(bams)}) =====")

            # --- samples found in the header
            try:
                samples = bam_samples(bam)
            except Exception as ex:
                self.log(f"[ERROR] Failed to read header: {ex}")
                fail.append(name)
                continue
            self.log(f"[INFO] Samples (SM) in the header: {', '.join(samples) if samples else 'none'}")

            tmp_bam = None
            if len(samples) == 1:
                jobs = [(bam, None, name)]
            elif len(samples) > 1 and cfg["multi"] == "split":
                jobs = [(bam, sm, f"{name}.{safe_name(sm)}") for sm in samples]
            else:
                # no SM tag, or merge requested
                if not samples:
                    self.log("[WARNING] No SM tag: read groups will be created with SM=" + name)
                out_final = os.path.join(cfg["out"], f"{name}.g.vcf.gz")
                if os.path.isfile(out_final) and not cfg["overwrite"]:
                    self.log(f"[INFO] Output already exists, skipped: {out_final}")
                    skip.append(name)
                    continue
                tmp_bam = self._unify_read_groups(cfg, bam, name)
                if not tmp_bam:
                    self.log("[ERROR] AddOrReplaceReadGroups failed.")
                    fail.append(name)
                    continue
                jobs = [(tmp_bam, None, name)]

            for in_bam, sm, label in jobs:
                if self.stop_flag.is_set():
                    break
                out_gvcf = os.path.join(cfg["out"], f"{label}.g.vcf.gz")
                logfile = os.path.join(cfg["out"], f"{label}.haplotypecaller.log")
                if os.path.isfile(out_gvcf) and not cfg["overwrite"]:
                    self.log(f"[INFO] Output already exists, skipped: {out_gvcf}")
                    skip.append(label)
                    continue
                if sm:
                    self.log(f"--- sample {sm}")
                ts = time.time()
                rc = self._exec(self._build_hc_cmd(cfg, in_bam, out_gvcf, sm), logfile)
                if self.stop_flag.is_set():
                    self.log("[INFO] Stopped by the user.")
                    for f in (out_gvcf, out_gvcf + ".tbi"):
                        if os.path.isfile(f):
                            os.remove(f)
                    break
                if rc == 0 and os.path.isfile(out_gvcf):
                    self.log(f"[OK] {out_gvcf} ({(time.time() - ts) / 60:.1f} min)")
                    ok.append(label)
                else:
                    self.log(f"[ERROR] HaplotypeCaller exited with code {rc}. Details: {logfile}")
                    fail.append(label)
            self._remove_tmp(tmp_bam)

        elapsed = (time.time() - t0) / 60
        self.log(f"\nSummary: {len(ok)} completed, {len(fail)} failed, {len(skip)} skipped — {elapsed:.1f} min")
        if fail:
            self.log("Failed: " + ", ".join(fail))
        self.msg_q.put(("progress", (len(bams), len(bams), "")))
        self.msg_q.put(("done", "Stopped" if self.stop_flag.is_set() else "Completed"))

    def _on_close(self):
        if self.proc and self.proc.poll() is None:
            if not messagebox.askyesno("Quit", "An analysis is in progress. Stop it and quit?"):
                return
            self._stop()
        self.destroy()


if __name__ == "__main__":
    GVCFApp().mainloop()
