# CYPredict

**CYPredict** is a machine learning-based web application that predicts whether drug-like molecules inhibit cytochrome P450 (CYP450) enzymes using machine learning models. It is designed to support students and researchers in understanding drug metabolism, and aims to reduce the need for animal testing in early drug development.

(images/header.png)


## Features
- Predict CYP inhibition from SMILES or CSV input
- Visualize molecule structure both 2D and 3D
- Display molecular properties
- View prediction confidence scores
- Download results in CSV/Excel


## Tech Stack
- Streamlit (UI)
- PaDEL-Descriptor (fingerprints)
- Scikit-learn (ML models)
- RDKit and PaDEL (molecule handling)



## Installation

### 1. Clone the Repository

### 2. Set Up the Environment

```bash
conda create --name cyp_pred python=3.8
conda activate cyp_pred
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```
### 4. Run the App
```bash
conda activate cyp_pred
streamlit run home.py
```






