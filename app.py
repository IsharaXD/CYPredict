import os
import subprocess
import streamlit as st
import pandas as pd
import joblib
import numpy as np
from functools import lru_cache

# ================== Setup ==================
st.set_page_config(
    page_title="Multi-CYP Inhibitor Predictor",
    page_icon="🧪",
    layout="wide"
)

# Custom CSS for three-column results
st.markdown("""
<style>
    .header {
        color: #4f8bf9;
        font-size: 24px !important;
    }
    .result-box {
        border-radius: 10px;
        padding: 15px;
        margin: 10px 0;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    .result-row {
        display: flex;
        gap: 10px;
    }
    .result-col {
        flex: 1;
        padding: 15px;
        border-radius: 8px;
    }
    .inhibitor {
        background-color: #853131;
    }
    .non-inhibitor {
        background-color: #3cc73c;
    }
    .probability-high {
        color: #ff0000;
        font-weight: bold;
    }
    .probability-medium {
        color: #ff9900;
    }
    .probability-low {
        color: #00aa00;
    }
</style>
""", unsafe_allow_html=True)

# ================== Descriptor Generator ==================
@lru_cache(maxsize=10)
def calculate_descriptors(smiles_text):
    try:
        smiles_list = smiles_text.strip().splitlines()
        with open('input.smi', 'w') as f:
            for smi in smiles_list:
                f.write(f"{smi.strip()}\t{smi.strip()}\n")

        cmd = [
            'java', '-Xms2G', '-Xmx2G', '-Djava.awt.headless=true',
            '-jar', 'PaDEL-Descriptor/PaDEL-Descriptor.jar',
            '-fingerprints', '-2d',
            '-descriptortypes', 'PaDEL-Descriptor/PubchemFingerprinter.xml',
            '-removesalt', '-standardizenitro', '-detectaromaticity',
            '-dir', './', '-file', 'descriptors_output.csv'
        ]

        result = subprocess.run(cmd, capture_output=True)
        if result.returncode != 0:
            st.error("⚠️ PaDEL descriptor generation failed.")
            st.code(result.stderr.decode())
            st.stop()

        desc = pd.read_csv('descriptors_output.csv')
        desc.columns = [col.strip() for col in desc.columns]
        return desc

    finally:
        for f in ['input.smi', 'descriptors_output.csv']:
            try: os.remove(f)
            except: pass

# ================== Model Loading ==================
@st.cache_resource
def load_models():
    try:
        # CYP1A2 Model
        cyp1a2_model = joblib.load('cyp1a2_xgboost.joblib')
        cyp1a2_features = joblib.load('padel_feature_names.pkl')
        
        # CYP2C19 Model
        cyp2c19_model = joblib.load('cyp2c19_xgboost.joblib')
        cyp2c19_features = joblib.load('padel_selected_features2c19.pkl')
        
        # CYP3A4 Model (NEW)
        cyp3a4_model = joblib.load('cyp3a4_xgboost.joblib')
        cyp3a4_features = joblib.load('padel_selected_features3a4.pkl')  # Ensure this has 267 features
        
        return {
            'cyp1a2': {
                'model': cyp1a2_model,
                'features': [str(f).strip() for f in cyp1a2_features],
                'labels': {0: "Inhibitor", 1: "Non-Inhibitor"}
            },
            'cyp2c19': {
                'model': cyp2c19_model,
                'features': [str(f).strip() for f in cyp2c19_features],
                'labels': {1: "Inhibitor", 0: "Non-Inhibitor"}
            },
            'cyp3a4': {  # NEW MODEL
                'model': cyp3a4_model,
                'features': [str(f).strip() for f in cyp3a4_features],
                'labels': {1: "Inhibitor", 0: "Non-Inhibitor"}  # Match your convention
            }
        }
    except Exception as e:
        st.error(f"❌ Model loading failed: {e}")
        st.stop()

# ================== UI Components ==================
def display_result(smiles, cyp1a2_pred, cyp2c19_pred, cyp3a4_pred):
    """Display results for one molecule with all three models"""
    # Determine CSS classes
    cyp1a2_class = "inhibitor" if cyp1a2_pred["prediction"] == "Inhibitor" else "non-inhibitor"
    cyp2c19_class = "inhibitor" if cyp2c19_pred["prediction"] == "Inhibitor" else "non-inhibitor"
    cyp3a4_class = "inhibitor" if cyp3a4_pred["prediction"] == "Inhibitor" else "non-inhibitor"
    
    # Probability styling
    def get_prob_class(prob):
        if prob > 0.7: return "probability-high"
        if prob > 0.5: return "probability-medium"
        return "probability-low"
    
    st.markdown(f"""
    <div class="result-box">
        <h3>SMILES: <code>{smiles}</code></h3>
        <div class="result-row">
            <div class="result-col {cyp1a2_class}">
                <h4>CYP1A2</h4>
                <p>Prediction: <strong>{cyp1a2_pred["prediction"]}</strong></p>
                <p>Probability: <span class="{get_prob_class(cyp1a2_pred["probability"])}">{cyp1a2_pred["probability"]:.3f}</span></p>
            </div>
            <div class="result-col {cyp2c19_class}">
                <h4>CYP2C19</h4>
                <p>Prediction: <strong>{cyp2c19_pred["prediction"]}</strong></p>
                <p>Probability: <span class="{get_prob_class(cyp2c19_pred["probability"])}">{cyp2c19_pred["probability"]:.3f}</span></p>
            </div>
            <div class="result-col {cyp3a4_class}">
                <h4>CYP3A4</h4>
                <p>Prediction: <strong>{cyp3a4_pred["prediction"]}</strong></p>
                <p>Probability: <span class="{get_prob_class(cyp3a4_pred["probability"])}">{cyp3a4_pred["probability"]:.3f}</span></p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

# ================== Main App ==================
st.title("🧪 Multi-CYP Inhibition Predictor")
st.markdown("""
Predict potential inhibition for **CYP1A2**, **CYP2C19**, and **CYP3A4** enzymes simultaneously.  
Enter SMILES strings below to analyze drug interaction risks.
""")

# Load models (cached)
models = load_models()

# Input section
with st.expander("🔬 Input SMILES", expanded=True):
    smiles_input = st.text_area(
        "Enter one or more SMILES strings (one per line):",
        height=150,
        value="CN1C=NC2=C1C(=O)N(C(=O)N2C)C\nCCOC(=O)c1ccccc1Cl\nCC(=O)OC1=CC=CC=C1C(=O)O"
    )
    predict_button = st.button("🚀 Predict for All", type="primary")

if predict_button and smiles_input:
    with st.spinner("🔍 Calculating predictions..."):
        try:
            # Generate descriptors
            desc = calculate_descriptors(smiles_input)
            smiles_list = smiles_input.strip().splitlines()
            
            results = []
            for i, smiles in enumerate(smiles_list):
                # Initialize feature matrices for all models
                cyp1a2_features = pd.DataFrame(0, index=[i], columns=models['cyp1a2']['features'])
                cyp2c19_features = pd.DataFrame(0, index=[i], columns=models['cyp2c19']['features'])
                cyp3a4_features = pd.DataFrame(0, index=[i], columns=models['cyp3a4']['features'])
                
                # Fill available features for each model
                for model_name in ['cyp1a2', 'cyp2c19', 'cyp3a4']:
                    features = models[model_name]['features']
                    available = [f for f in features if f in desc.columns]
                    
                    if model_name == 'cyp1a2':
                        cyp1a2_features[available] = desc.iloc[i][available]
                    elif model_name == 'cyp2c19':
                        cyp2c19_features[available] = desc.iloc[i][available]
                    else:  # cyp3a4
                        cyp3a4_features[available] = desc.iloc[i][available]
                
                # Make predictions
                cyp1a2_probs = models['cyp1a2']['model'].predict_proba(cyp1a2_features)[0]
                cyp2c19_probs = models['cyp2c19']['model'].predict_proba(cyp2c19_features)[0]
                cyp3a4_probs = models['cyp3a4']['model'].predict_proba(cyp3a4_features)[0]
                
                # Store results
                results.append({
                    "smiles": smiles,
                    "cyp1a2": {
                        "probability": cyp1a2_probs[0],
                        "prediction": models['cyp1a2']['labels'][np.argmax(cyp1a2_probs)]
                    },
                    "cyp2c19": {
                        "probability": cyp2c19_probs[1],
                        "prediction": models['cyp2c19']['labels'][np.argmax(cyp2c19_probs)]
                    },
                    "cyp3a4": {
                        "probability": cyp3a4_probs[1],
                        "prediction": models['cyp3a4']['labels'][np.argmax(cyp3a4_probs)]
                    }
                })
            
            # Display results
            st.success(f"✅ Successfully predicted {len(results)} molecules")
            
            # Summary table
            st.subheader("📊 Prediction Summary")
            summary_data = []
            for res in results:
                summary_data.append({
                    "SMILES": res["smiles"],
                    "CYP1A2 Prediction": res["cyp1a2"]["prediction"],
                    "CYP1A2 Probability": f"{res['cyp1a2']['probability']:.3f}",
                    "CYP2C19 Prediction": res["cyp2c19"]["prediction"],
                    "CYP2C19 Probability": f"{res['cyp2c19']['probability']:.3f}",
                    "CYP3A4 Prediction": res["cyp3a4"]["prediction"],
                    "CYP3A4 Probability": f"{res['cyp3a4']['probability']:.3f}"
                })
            
            st.dataframe(pd.DataFrame(summary_data), use_container_width=True)
            
            # Detailed results
            st.subheader("🔍 Detailed Results")
            for res in results:
                display_result(res["smiles"], res["cyp1a2"], res["cyp2c19"], res["cyp3a4"])
            
            # Download option
            csv = pd.DataFrame(summary_data).to_csv(index=False)
            st.download_button(
                "💾 Download All Results",
                data=csv,
                file_name="cyp_predictions.csv",
                mime="text/csv"
            )
            
        except Exception as e:
            st.error(f"❌ Prediction failed: {str(e)}")

# Add info sidebar
with st.sidebar:
    st.header("About This Tool")
    st.markdown("""
    This tool predicts inhibition potential for:
    - **CYP1A2**: {f1} features
    - **CYP2C19**: {f2} features
    - **CYP3A4**: {f3} features (metabolizes ~50% of drugs)
    
    **Interpretation Guide**:
    - 🟢 **Non-Inhibitor** (Probability < 0.5)
    - 🟠 **Potential Inhibitor** (0.5 ≤ Probability < 0.7)
    - 🔴 **Likely Inhibitor** (Probability ≥ 0.7)
    """.format(
        f1=len(models['cyp1a2']['features']),
        f2=len(models['cyp2c19']['features']),
        f3=len(models['cyp3a4']['features'])
    ))
    
    st.header("Example SMILES")
    st.code("Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    st.code("Omeprazole: CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=CC(=C3)OC")
    st.code("Midazolam: CN1C(=O)CN=C(C2=C1C=CC(=N2)C1=CC=CC=C1)C1=CC=CC=C1F")