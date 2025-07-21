# MoleculeAI 🔬💊

## 🧬 Introduction

MoleculeAI is an intelligent platform designed to assist in **drug discovery** using advanced machine learning models. It generates novel molecules, retrieves bioactivity data, and predicts their drug-likeness through ADMET analysis, ultimately helping researchers find potential candidates faster and more efficiently.

---

## ⚙️ Tech Stack

- **Frontend:** React.js, TailwindCSS
- **Backend:** Python (Flask / FastAPI)
- **ML Libraries:** scikit-learn, RDKit, TensorFlow/PyTorch (VAE, Fragment-Based Generator)
- **Data Sources:** ChEMBL, STRING DB
- **Deployment:** Vercel

---

## 🚀 Features

- 🔍 Search compounds by protein targets (e.g., DRD2, EGFR, BRAF)
- 📡 Automatic bioactivity data retrieval from ChEMBL
- 🧹 Preprocessing pipeline for cleaning and fingerprint generation
- 🧠 AI-powered molecule generation (VAE + Fragment Generator)
- 💊 ADMET property filtering
- 📋 Summary report generation
- 🌙 Light/Dark Mode toggle UI
- 💻 Clean and responsive UI for researchers and developers

---

## ⚡ Quick Start

### 🔧 Prerequisites

- Node.js and npm installed
- Python 3.8+
- Git

---

### 📥 Clone the Repository

```bash
git clone https://github.com/yourusername/moleculeai.git
cd moleculeai

```
```bash
cd backend
pip install -r requirements.txt
python app.py
```

```bash
cd frontend
npm install
npm run dev
```

