# Gadgetron Parallel Reconstruction Framework

**A practical framework for translating offline MR reconstructions to inline deployment, built on the open-source Gadgetron platform.**

---

## 🧭 Overview

This repository provides a plug-and-play framework for converting custom MR reconstructions—originally developed for offline use—into robust, scalable, and clinically usable inline workflows.

* ✔️ Compatible with **Siemens platforms**
* 🔁 Enables **multi-input**, **parallel & resource-aware**, and **asynchronous** inline execution
* 🧩 Preserves **scanner-based post-processing** for unified image appearance and clinical review

This framework has been validated in both routine clinical protocols and large-scale studies such as **TwinsUK**, supporting inline deployment of reconstructions like **SENSE**, **AlignSENSE**, and **radial NUFFT**.

---

## ✨ Key Features

### ✅ Easy Prototyping of Offline Reconstructions

* Converts ISMRMRD raw data into a **Twix-like structure** for minimal modification of existing offline code
* Provides a general-purpose **MATLAB-based input converter** and customizable **ParameterMap** covering commonly used headers

### ⚙️ Inline Execution Without Scanner Disruption

* Enables **long custom reconstructions** without interrupting the MR exam flow
* Inline reconstructions are executed **asynchronously** on an external server
* Supports **retrieval scans** and **retro-reconstruction**
* Scanner-native reconstructions are preserved for **early review**

### 🔁 Multi-input & Joint Reconstruction Support

* Allows inline reconstructions to access **shared or external input scans** (e.g., reference coils)

### 🧠 Resource-aware Parallel Scheduling

* Includes built-in GPU monitoring and queueing
* Schedules jobs intelligently across **multi-GPU servers** to avoid overloads

### 🧬 Integrated Scanner Post-processing

* Enables scanner-based post-processing (e.g., bias field and distortion correction) on custom reconstructions
* Ensures **visual consistency** with Siemens-reconstructed images

---

## 📁 Folder Structure

```
Gadgetron_Parallel_Framework/
│
├── ParameterMap/             # Metadata mapping (Twix → ISMRMRD)
├── config/                   # Gadget chain configuration files (*.xml)
├── matlab_handle_archive/    # MATLAB handlers for "Read & Save" and background recon
├── Gadgetron_tools/          # Input converters, Twix-like wrappers, utilities
├── Bash_script/              # Headless / queued job launchers
├── Useful_tools/             # Misc tools: ISMRMRD readers, resource monitors, etc.
├── SENSE_Recon_Demo/         # Demo: SENSE (offline + inline)
├── AlignSENSE_Recon_Demo/    # Demo: AlignSENSE (offline + inline)
```

---

## 🚀 Getting Started

### Requirements

* [Gadgetron](https://github.com/gadgetron/gadgetron)
* **MATLAB ≥ R2021a** with Gadgetron Toolbox
* **CUDA-enabled Conda environment** (for launching recon scripts)
* **Linux server with NVIDIA GPUs**
* **Siemens scanner with ICEGadgetron** integration

### 📖 Full Documentation

All setup instructions and tutorials are available at:

👉 [User Manual on Notion](https://shine-pond-caf.notion.site/User-manual-20961ff38021807f89a8fdcc819acd0b?source=copy_link)

---

## 🧪 Demo Cases

Inline implementations included:

* **SENSE**: `SENSE_Recon_Demo/`
* **AlignSENSE**: `AlignSENSE_Recon_Demo/`
* **radial NUFFT** (under `matlab_handle_archive/`)

Each demo contains:

* Offline version
* Inline version (`inline/` subfolder)
* `diff/` subfolder showing minimal modifications for inline adaptation

---

## ✅ Validation Summary

Validated on:

* **3T MAGNETOM Vida (XA60)** — tested with MPRAGE, FLAIR, SWI
* **7T MAGNETOM Terra X** — tested with multi-echo radial GRE
* **TwinsUK cohort study**

  * N = 480 subjects
  * 99% successful inline retrieval rate (via retrieval scans or retro-recon)

---

## ⚠️ Limitations and Roadmap

* ✅ Siemens only (FIRE/OpenRecon compatibility under discussion)
* 🧪 **MATLAB-only input converter** (Python version coming soon)
* 📦 Raw data currently saved as `.mat` (future HDF5 support planned)

---

## 📝 Citation

If you use this framework in your research, please cite:

> Ning Z, et al. *"From Offline to Inline Without Pain: A Practical Framework for Translating Offline MR Reconstructions to Inline Deployment Using the Gadgetron Platform."* (2025)

---

## 📄 License

This project is released under the [MIT License](LICENSE).
