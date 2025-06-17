# Gadgetron Parallel Reconstruction Framework

A MATLAB-based extension to the Gadgetron medical-image reconstruction platform that enables headless, parallelizable offline and inline workflows.

## Features

- **Offline Test Modes**  
  - **Debug Mode**: Launch MATLAB with full GUI, set breakpoints, and step through data handlers.  
  - **Auto Mode**: Headless execution via the Gadgetron client, with resource-aware handoff to background reconstruction jobs.

- **Inline Deployment**  
  - Integrates with Siemens ICE (XA30, XA60, VE12) via IceGadgetron.  
  - Configurable ParameterMaps, ICE workflows, and TCP/IP connectivity to your external server.

- **Built-in Demos**  
  - SENSE, AlignSENSE, and radial NUFFT reconstructions.  
  - MATLAB handler templates for “Read & Save” and “Background Recon” stages.  
  - Automated GPU-resource monitoring and queuing via Bash scripts.

## Prerequisites

See the [User Manual](https://shine-pond-caf.notion.site/User-manual-20961ff38021807f89a8fdcc819acd0b?pvs=73) for full details, but at minimum you’ll need:

- **Gadgetron client** (Conda environment with CUDA)  
- **MATLAB** (R2021a or later) with the Gadgetron MATLAB Toolbox  
- **Siemens ICE** environment and IceGadgetron binaries (requires a Siemens–Gadgetron agreement)  
- Network access between the MR console (MARS) and your reconstruction server

## Quick Start

1. **Clone the repository**  
```
git clone https://github.com/ZihanNing/Gadgetron_Parallel_Framework.git
cd Gadgetron_Parallel_Framework
````

2. **Refer to the User Manual**
   Follow the step-by-step instructions for:

   * Offline Debug Mode
   * Offline Auto Mode
   * Inline Deployment

   👉 [Open the User Manual](https://shine-pond-caf.notion.site/User-manual-20961ff38021807f89a8fdcc819acd0b?pvs=73)

---

© 2025 ZihanNing | MIT License
