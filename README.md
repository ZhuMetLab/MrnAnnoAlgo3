# MrnAnnoAlgo3 (MetDNA3) <img src="man/figures/logo.png" align="right" alt="MetDNA3 Logo" width="150"/>

[![GitHub release](https://img.shields.io/github/v/release/ptpb781209069/MrnAnnoAlgo3)](https://github.com/ptpb781209069/MrnAnnoAlgo3)
[![License](https://img.shields.io/badge/license-GPL--3-blue)](https://opensource.org/licenses/GPL-3.0)

**Knowledge- and Data-Driven Metabolite Annotation Engine for Metabolomics**

`MrnAnnoAlgo3` is the core algorithm module of **MetDNA3**, designed to annotate metabolites through a two-layer interactive networking topology (knowledge- and data-driven) and recursive annotation algorithms. It provides a robust computational foundation for large-scale metabolomic studies.

---

## Table of Contents
- [Key Features](#key-features)
- [Installation](#installation)
- [Usage](#usage)
- [Documentation](#documentation)
- [Support](#support)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

---

## Key Features
✨ **Dual-Layer Annotation**  
Combines **knowledge-driven** (e.g., biochemical pathways, spectral libraries) and **data-driven** (experimental MS/MS networks) layers for comprehensive metabolite annotation.

⚡ **Recursive Algorithm**  
Iteratively refines annotation confidence by propagating spectral similarity across the network.

📊 **High Performance**  
Optimized for large-scale datasets with efficient C++ backend (via `Rcpp`).

🔗 **Seamless Integration**  
Part of the MetDNA3 ecosystem. *Note: Full MetDNA3 functionality requires additional modules.*

---

## Installation

### From GitHub (Recommended)
```r
# Install via devtools (ensure devtools is installed)
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ptpb781209069/MrnAnnoAlgo3")
```

### From Gitea (Internal/Mirror)
```r
devtools::install_git("http://10.42.0.33:3000/zhanghs/MrnAnnoAlgo3.git")
```

---

## Usage
### Basic Example
```r
library(MrnAnnoAlgo3)

# Load example data
data(sample_spectra)

# Initialize annotation engine
annotator <- new("MetDNA3Annotator", 
                 knowledge_network = "path_to_knowledge_graph",
                 experimental_data = sample_spectra)

# Run recursive annotation
results <- annotateMetabolites(annotator)
```

### Advanced Configuration
See the [Documentation Website](https://metdna3-docs.example.com) for detailed tutorials on:
- Customizing knowledge networks
- Tuning recursive algorithm parameters
- Interpreting annotation confidence scores

---

## Documentation
- 📚 **Full Documentation**: [MetDNA3 Algorithm Guide](https://metdna3-docs.example.com/mrnannoalgo3)
- 📄 **Vignettes**: Run `browseVignettes("MrnAnnoAlgo3")` after installation.

---

## Support
- 🐛 **Bug Reports & Feature Requests**: [GitHub Issues](https://github.com/ptpb781209069/MrnAnnoAlgo3/issues)
- 📧 **Direct Contact**: zhanghs@sioc.ac.cn
- 💬 **Community Forum**: [MetDNA3 Discussions](https://github.com/orgs/MetDNA3/discussions) (Coming Soon)

---

## Citation
If you use `MrnAnnoAlgo3` in your research, please cite:

> Zhang H, et al. (2024). *MetDNA3: A Scalable Framework for Knowledge- and Data-Driven Metabolite Annotation*. *Nature Methods* (Under Review). DOI: [Pending]

*Check back later for updates or email for preprint access.*

---

## Contributing
We welcome contributions! Please follow these steps:
1. Fork the repository.
2. Create a branch (`git checkout -b feature/your-feature`).
3. Commit changes (`git commit -am 'Add some feature'`).
4. Push to the branch (`git push origin feature/your-feature`).
5. Open a [Pull Request](https://github.com/ptpb781209069/MrnAnnoAlgo3/pulls).

---

## License
This project is licensed under the **GNU General Public License v3.0**. See [LICENSE](LICENSE) for details.
