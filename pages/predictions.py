import os
import subprocess
import streamlit as st
import pandas as pd
import joblib
import numpy as np
from functools import lru_cache
from app import show_navbar
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
from io import BytesIO
import cairosvg
import requests
import py3Dmol
import base64
from io import BytesIO

# ================== Setup ==================
st.set_page_config(
    page_title="Multi-CYP Inhibitor Predictor",
    page_icon="🧪",
    layout="wide"
)

selected = show_navbar()

# Custom CSS for improved structure
st.markdown("""
<style>
    /* Main containers */
    .main-container {
        padding: 20px;
        border-radius: 10px;
        background-color: #f9f9f9;
        margin-bottom: 20px;
    }
    .molecule-container {
        padding: 15px;
        border-radius: 8px;
        background-color: white;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        margin-bottom: 15px;
    }
    
    /* Typography */
    .section-title {
        color: #2c3e50;
        border-bottom: 2px solid #3498db;
        padding-bottom: 5px;
        margin-bottom: 15px;
    }
    .subsection-title {
        color: #3498db;
        margin-top: 15px;
        margin-bottom: 10px;
    }
    
    /* Results styling */
    .result-card {
        border-radius: 8px;
        padding: 15px;
        margin: 10px 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        background-color: white;
    }
    .result-row {
        display: flex;
        gap: 15px;
        margin-bottom: 15px;
    }
    .result-col {
        flex: 1;
        padding: 15px;
        border-radius: 8px;
        text-align: center;
    }
    .inhibitor {
        background-color: #ffebee;
        border-left: 4px solid #f44336;
    }
    .non-inhibitor {
        background-color: #e8f5e9;
        border-left: 4px solid #4caf50;
    }
    .probability-high {
        color: #f44336;
        font-weight: bold;
    }
    .probability-medium {
        color: #ff9800;
    }
    .probability-low {
        color: #4caf50;
    }
    
    /* Tables */
    .property-table {
        width: 100%;
        border-collapse: collapse;
    }
    .property-table th, .property-table td {
        padding: 10px;
        border: 1px solid #ddd;
        text-align: left;
    }
    .property-table th {
        background-color: #f2f2f2;
    }
    
    /* Navigation */
    .tab-container {
        margin-bottom: 20px;
    }
</style>
""", unsafe_allow_html=True)

# ================== Helper Functions ==================
def validate_smiles(smiles):
    """Validate SMILES string"""
    mol = Chem.MolFromSmiles(smiles)
    return (True, mol) if mol else (False, "Invalid SMILES string")

def get_compound_names(smiles):
    """Fetch IUPAC and common names"""
    try:
        iupac = requests.get(f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/iupac_name", timeout=3).text.strip()
    except:
        iupac = None
    
    common_name, synonyms = None, []
    try:
        response = requests.get(f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/names", timeout=3)
        if response.ok:
            names = [n.strip() for n in response.text.split('\n') if n.strip()]
            filtered = [n for n in names if not any(c.isdigit() for c in n[:4])]
            common_name = filtered[0] if filtered else None
            synonyms = filtered or names
    except:
        pass
    
    return {
        "iupac": iupac or "Not available",
        "common_name": common_name or "Unnamed compound",
        "synonyms": synonyms[:10]  # Limit to 10 synonyms
    }

def draw_2d_molecule(mol, size=(400, 400)):
    """Generate 2D molecule image"""
    try:
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(*size)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        png = cairosvg.svg2png(bytestring=svg.encode())
        return Image.open(BytesIO(png))
    except Exception as e:
        st.error(f"2D rendering error: {e}")
        return None

def show_3d_molecule(smiles):
    """Generate interactive 3D visualization"""
    try:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        viewer = py3Dmol.view(width=600, height=400)
        viewer.addModel(Chem.MolToMolBlock(mol), 'mol')
        viewer.setStyle({'stick': {}})
        viewer.zoomTo()
        return viewer
    except Exception as e:
        st.error(f"3D rendering error: {e}")
        return None

def get_molecule_properties(mol):
    """Calculate molecular properties"""
    return {
        "Molecular Weight": f"{Descriptors.MolWt(mol):.2f} g/mol",
        "LogP": f"{Descriptors.MolLogP(mol):.2f}",
        "H-Bond Donors": Descriptors.NumHDonors(mol),
        "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": f"{Descriptors.TPSA(mol):.2f} Å²",
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Aromatic Rings": Descriptors.NumAromaticRings(mol),
        "Heavy Atoms": Descriptors.HeavyAtomCount(mol)
    }

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
        cyp1a2_model = joblib.load('models/cyp1a2_xgboost.joblib')
        cyp1a2_features = joblib.load('pkl/padel_feature_names.pkl')
        
        # CYP2C19 Model
        cyp2c19_model = joblib.load('models/cyp2c19_xgboost.joblib')
        cyp2c19_features = joblib.load('pkl/padel_selected_features2c19.pkl')
        
        # CYP3A4 Model
        cyp3a4_model = joblib.load('models/cyp3a4_xgboost.joblib')
        cyp3a4_features = joblib.load('pkl/padel_selected_features3a4.pkl')
        
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
            'cyp3a4': {
                'model': cyp3a4_model,
                'features': [str(f).strip() for f in cyp3a4_features],
                'labels': {1: "Inhibitor", 0: "Non-Inhibitor"}
            }
        }
    except Exception as e:
        st.error(f"❌ Model loading failed: {e}")
        st.stop()

# ================== UI Components ==================
def display_molecule_card(smiles):
    """Display molecular information in a structured card"""
    is_valid, mol = validate_smiles(smiles)
    if not is_valid:
        st.error(f"Invalid SMILES: {smiles}")
        return None
    
    names = get_compound_names(smiles)
    props = get_molecule_properties(mol)
    img_2d = draw_2d_molecule(mol)
    
    with st.container():
        st.markdown(f"### {names['common_name']}")
        st.caption(f"IUPAC: {names['iupac']}")
        
        # Main columns
        col1, col2 = st.columns([1, 2])
        
        with col1:
            # 2D Structure
            if img_2d:
                st.image(img_2d, use_container_width=True)
            
            # Properties table
            st.markdown("**Molecular Properties**")
            props_df = pd.DataFrame(list(props.items()), columns=["Property", "Value"])
            st.dataframe(props_df, hide_index=True, use_container_width=True)
        
        with col2:
            # 3D Structure
            st.markdown("**3D Visualization**")
            viewer = show_3d_molecule(smiles)
            if viewer:
                st.components.v1.html(viewer._make_html(), height=400)
            
            # Synonyms
            if names['synonyms']:
                with st.expander(f"Synonyms ({len(names['synonyms'])})"):
                    for syn in names['synonyms']:
                        st.markdown(f"- {syn}")
    
    return mol

def display_prediction_result(smiles, preds):
    """Display prediction results in a structured card"""
    # Determine CSS classes
    cyp1a2_class = "inhibitor" if preds["cyp1a2"]["prediction"] == "Inhibitor" else "non-inhibitor"
    cyp2c19_class = "inhibitor" if preds["cyp2c19"]["prediction"] == "Inhibitor" else "non-inhibitor"
    cyp3a4_class = "inhibitor" if preds["cyp3a4"]["prediction"] == "Inhibitor" else "non-inhibitor"
    
    # Probability styling
    def get_prob_class(prob):
        if prob > 0.7: return "probability-high"
        if prob > 0.5: return "probability-medium"
        return "probability-low"
    
    with st.container():
        st.markdown(f"### Prediction Results for `{smiles}`")
        
        # Results in columns
        cols = st.columns(3)
        
        with cols[0]:
            st.markdown(f"<div class='result-col {cyp1a2_class}'>"
                       f"<h4>CYP1A2</h4>"
                       f"<p><strong>{preds['cyp1a2']['prediction']}</strong></p>"
                       f"<p class='{get_prob_class(preds['cyp1a2']['probability'])}'>"
                       f"Probability: {preds['cyp1a2']['probability']:.3f}</p>"
                       "</div>", unsafe_allow_html=True)
        
        with cols[1]:
            st.markdown(f"<div class='result-col {cyp2c19_class}'>"
                       f"<h4>CYP2C19</h4>"
                       f"<p><strong>{preds['cyp2c19']['prediction']}</strong></p>"
                       f"<p class='{get_prob_class(preds['cyp2c19']['probability'])}'>"
                       f"Probability: {preds['cyp2c19']['probability']:.3f}</p>"
                       "</div>", unsafe_allow_html=True)
        
        with cols[2]:
            st.markdown(f"<div class='result-col {cyp3a4_class}'>"
                       f"<h4>CYP3A4</h4>"
                       f"<p><strong>{preds['cyp3a4']['prediction']}</strong></p>"
                       f"<p class='{get_prob_class(preds['cyp3a4']['probability'])}'>"
                       f"Probability: {preds['cyp3a4']['probability']:.3f}</p>"
                       "</div>", unsafe_allow_html=True)

def create_comprehensive_report(results, smiles_list):
    report_data = []
    
    for i, (smiles, res) in enumerate(zip(smiles_list, results)):
        # Get molecular information
        is_valid, mol = validate_smiles(smiles)
        if not is_valid:
            continue
            
        names = get_compound_names(smiles)
        props = get_molecule_properties(mol)
        
        # Create 2D structure image (base64 encoded)
        img_2d = None
        try:
            img = draw_2d_molecule(mol, size=(300, 300))
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_2d = base64.b64encode(buffered.getvalue()).decode()
        except:
            img_2d = ""
        
        # Compile all data
        report_data.append({
            "Molecule_ID": f"Molecule_{i+1}",
            "SMILES": smiles,
            "Common_Name": names['common_name'],
            "IUPAC_Name": names['iupac'],
            "Synonyms": "; ".join(names['synonyms']),
            **props,  # Unpacks all molecular properties
            "CYP1A2_Prediction": res["cyp1a2"]["prediction"],
            "CYP1A2_Probability": res["cyp1a2"]["probability"],
            "CYP2C19_Prediction": res["cyp2c19"]["prediction"],
            "CYP2C19_Probability": res["cyp2c19"]["probability"],
            "CYP3A4_Prediction": res["cyp3a4"]["prediction"],
            "CYP3A4_Probability": res["cyp3a4"]["probability"],
            "2D_Structure_Base64": img_2d
        })
    
    return pd.DataFrame(report_data)

# ================== Main App ==================
def main():
    st.title("🧪 Multi-CYP Inhibition Predictor")
    st.markdown("""
    Predict potential inhibition for **CYP1A2**, **CYP2C19**, and **CYP3A4** enzymes.  
    Enter SMILES strings below to analyze drug interaction risks.
    """)
    
    # Load models
    models = load_models()
    
    # Input section
    with st.expander("🔬 Input SMILES", expanded=True):
        smiles_input = st.text_area(
            "Enter one or more SMILES strings (one per line):",
            height=150,
            value="CN1C=NC2=C1C(=O)N(C(=O)N2C)C\nCCOC(=O)c1ccccc1Cl\nCC(=O)OC1=CC=CC=C1C(=O)O",
            help="Example: Caffeine (CN1C=NC2=C1C(=O)N(C(=O)N2C)C)"
        )
        predict_button = st.button("🚀 Predict", type="primary")
    
    if predict_button and smiles_input:
        with st.spinner("🔍 Calculating predictions..."):
            try:
                # Generate descriptors
                desc = calculate_descriptors(smiles_input)
                smiles_list = smiles_input.strip().splitlines()
                
                results = []
                for i, smiles in enumerate(smiles_list):
                    # Display molecule information
                    with st.container():
                        st.markdown(f"## Molecule {i+1}")
                        mol = display_molecule_card(smiles)
                        if not mol:
                            continue
                    
                    # Initialize feature matrices
                    cyp1a2_features = pd.DataFrame(0, index=[i], columns=models['cyp1a2']['features'])
                    cyp2c19_features = pd.DataFrame(0, index=[i], columns=models['cyp2c19']['features'])
                    cyp3a4_features = pd.DataFrame(0, index=[i], columns=models['cyp3a4']['features'])
                    
                    # Fill available features
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
                    preds = {
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
                    }
                    
                    # Display prediction results
                    display_prediction_result(smiles, preds)
                    results.append({"smiles": smiles, **preds})
                    
                    # Divider between molecules
                    if i < len(smiles_list) - 1:
                        st.divider()
                
                # Summary and download
                if results:
                    st.markdown("<br>", unsafe_allow_html=True)  # Add space above the box
                    st.success(f"✅ Successfully predicted {len(results)} molecules")
                    # Create two tabs for different views
                    tab1, tab2 = st.tabs(["📊 Summary Table", "💾 Download Options"])
                    
                    with tab1:
                        # Summary table view
                        summary_data = []
                        for res in results:
                            summary_data.append({
                                "SMILES": res["smiles"],
                                "CYP1A2": f"{res['cyp1a2']['prediction']} ({res['cyp1a2']['probability']:.3f})",
                                "CYP2C19": f"{res['cyp2c19']['prediction']} ({res['cyp2c19']['probability']:.3f})",
                                "CYP3A4": f"{res['cyp3a4']['prediction']} ({res['cyp3a4']['probability']:.3f})"
                            })
                        st.dataframe(pd.DataFrame(summary_data), use_container_width=True)
                    
                    with tab2:
                        # Download options
                        st.subheader("Download Options")
                        
                        # Option 1: Simple CSV (just predictions)
                        st.markdown("**Basic Prediction Data (CSV)**")
                        csv = pd.DataFrame(summary_data).to_csv(index=False)
                        st.download_button(
                            "⬇️ Download Predictions (CSV)",
                            data=csv,
                            file_name="cyp_predictions.csv",
                            mime="text/csv",
                            help="Simple table with just prediction results"
                        )
                        
                        # Option 2: Full Excel Report
                        st.markdown("**Comprehensive Report (Excel)**")
                        with st.spinner("Preparing full report..."):
                            report_df = create_comprehensive_report(results, smiles_list)
                            
                            output = BytesIO()
                            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                                # Main data sheet
                                report_df.drop(columns=['2D_Structure_Base64']).to_excel(writer, sheet_name='Results', index=False)
                                
                                # Structures sheet
                                workbook = writer.book
                                worksheet = workbook.add_worksheet('Structures')
                                
                                # Add images to Excel
                                for i, img_data in enumerate(report_df['2D_Structure_Base64']):
                                    if img_data:
                                        worksheet.insert_image(
                                            f'A{i*15+1}',  # Add spacing between images
                                            f'Structure_{i+1}.png', 
                                            {'image_data': BytesIO(base64.b64decode(img_data))}
                                        )
                                
                                # Add metadata sheet
                                metadata = pd.DataFrame({
                                    'Column': report_df.columns,
                                    'Description': [
                                        'Sequential molecule ID',
                                        'SMILES string',
                                        'Common chemical name',
                                        'IUPAC systematic name',
                                        'Semicolon-separated list of synonyms',
                                        *[f'{k} (molecular property)' for k in get_molecule_properties(Chem.MolFromSmiles("C")).keys()],
                                        'CYP1A2 inhibition prediction',
                                        'CYP1A2 prediction probability (0-1)',
                                        'CYP2C19 inhibition prediction',
                                        'CYP2C19 prediction probability (0-1)',
                                        'CYP3A4 inhibition prediction',
                                        'CYP3A4 prediction probability (0-1)',
                                        'Base64 encoded PNG structure (internal use)'
                                    ]
                                })
                                metadata.to_excel(writer, sheet_name='Metadata', index=False)
                            
                            excel_data = output.getvalue()
                        
                        st.download_button(
                            "⬇️ Download Full Report (Excel)",
                            data=excel_data,
                            file_name="cyp_predictions_full_report.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            help="Includes all molecular information, properties, and structures"
                        )
            
            except Exception as e:
                st.error(f"❌ Prediction failed: {str(e)}")

# Sidebar
st.markdown("""
    <style>
        section[data-testid="stSidebar"] {
            display: block !important;
            width: 350px !important;
        }
    </style>
    """, unsafe_allow_html=True)
    
    # 3. Now create your sidebar content
with st.sidebar:
        # Hide auto-generated pages list (but keep sidebar visible)
        st.markdown("""
        <style>
        /* Show sidebar */
         section[data-testid="stSidebar"] {
            display: block !important;
            width: 350px !important;
        }
        
        /* Hide auto-generated pages list */
        div[data-testid="stSidebarNav"] {
            display: none !important;
    </style>
    """, unsafe_allow_html=True)
    
    # 3. Your custom sidebar content
with st.sidebar:
        st.header("ℹ️ About This Tool")
        models = load_models()
        
        st.markdown(f"""
        Predicts inhibition potential for key CYP enzymes:
        - **CYP1A2**: {len(models['cyp1a2']['features'])} features
        - **CYP2C19**: {len(models['cyp2c19']['features'])} features
        - **CYP3A4**: {len(models['cyp3a4']['features'])} features
        """)
        
        st.markdown("""
        **Interpretation Guide**:
        - 🟢 **Non-Inhibitor** (Probability < 0.5)
        - 🟠 **Potential Inhibitor** (0.5 ≤ Probability < 0.7)
        - 🔴 **Likely Inhibitor** (Probability ≥ 0.7)
        """)
        
        st.header("🧪 Example SMILES")
        st.code("Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        st.code("Omeprazole: CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=CC(=C3)OC")
        st.code("Midazolam: CN1C(=O)CN=C(C2=C1C=CC(=N2)C1=CC=CC=C1)C1=CC=CC=C1F")

    # Main content area
    # ... your prediction content ...

if __name__ == "__main__":
    main()