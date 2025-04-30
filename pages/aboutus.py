import streamlit as st
from PIL import Image
import os
from app import show_navbar

# Custom CSS with expanded layout
st.set_page_config(layout="wide")
st.markdown("""
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
<style>
    .main {
        max-width: 1200px;
        padding: 2rem 4rem;
    }
    .section {
        padding: 2rem;
        border-radius: 15px;
        margin: 1.5rem 0;
        background-color: #f8f9fa;
        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        transition: all 0.3s ease;
    }
    .section:hover {
        box-shadow: 0 8px 16px rgba(0,0,0,0.15);
    }
    .profile-img-container {
        display: flex;
        justify-content: center;
        margin-bottom: 1.5rem;
    }
    .profile-img {
        border-radius: 50%;
        border: 4px solid #4a8cff;
        object-fit: cover;
        width: 250px;
        height: 250px;
    }
    .social-icon {
        font-size: 2rem;
        color: #4a8cff;
        margin-bottom: 0.5rem;
    }
    h2 {
        color: #2c3e50;
        border-bottom: 2px solid #4a8cff;
        padding-bottom: 0.5rem;
        margin-top: 0;
    }
    .tech-card {
        background: white;
        border-radius: 10px;
        padding: 1.5rem;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        height: 100%;
    }
    .contact-card {
        text-align: center;
        padding: 1.5rem;
        border-radius: 10px;
        transition: all 0.3s ease;
        background: white;
    }
    .contact-card:hover {
        transform: translateY(-5px);
        box-shadow: 0 6px 12px rgba(0,0,0,0.15);
    }
    .grid-container {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
        gap: 1.5rem;
    }
    @media (max-width: 768px) {
        .main {
            padding: 1rem;
        }
    }
</style>
""", unsafe_allow_html=True)

selected = show_navbar()

def about_page():
    # Load profile image
    try:
        profile_img = Image.open("images/me.jpeg")
    except:
        profile_img = Image.new('RGB', (250, 250), color='lightgray')
        st.warning("Profile image not found - using placeholder")

    
    
    # Hero Section
    col1, col2 = st.columns([1, 2], gap="large")
    
    with col1:
        st.markdown("""
        <div class="profile-img-container" style="margin-top: 20px;">
            <img src="data:image/png;base64,{}" class="profile-img" style="border: none; box-shadow: 0 4px 8px rgba(0,0,0,0.2);">
        </div>
        """.format(image_to_base64(profile_img)), unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class='section'>
        <h2><i class="fas fa-user" ></i> About Me</h2>
        <p style='font-size:1.1rem; line-height:1.6;'>
        Hi! I'm Pallawattha Bhagya, a final-year BSc (Hons) Software Engineering student at NSBM Green University. 
        I am passionate about ethical, science-driven technologies, especially those that make education and 
        research more accessible.
        </p>
        <p style='font-size:1.1rem; line-height:1.6;'>
        This project, "CYP450 Inhibition Prediction Tool", was developed as part of my final year individual 
        research project. It combines my interests in machine learning, cheminformatics, and ethical drug 
        development to build a user-friendly tool that can predict enzyme inhibition.
        </p>
        </div>
        """, unsafe_allow_html=True)
    
    # What I Do Section
    st.markdown("""
    <div class='section'>
    <h2><i class="fas fa-brain"></i> What I Do</h2>
    <div class="grid-container">
    <div class="tech-card">
    <h4><i class="fas fa-robot"></i> AI Development</h4>
    <p>Develop ethical AI tools for scientific research and drug discovery</p>
    </div>
    <div class="tech-card">
    <h4><i class="fas fa-atom"></i> Cheminformatics</h4>
    <p>Work with SMILES, molecular fingerprints and chemical data</p>
    </div>
    <div class="tech-card">
    <h4><i class="fas fa-laptop-code"></i> UI Development</h4>
    <p>Build clean, intuitive interfaces with Streamlit</p>
    </div>
    <div class="tech-card">
    <h4><i class="fas fa-lightbulb"></i> Problem Solving</h4>
    <p>Translate scientific challenges into software solutions</p>
    </div>
    </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Mission Section
    st.markdown("""
    <div class='section'>
    <h2><i class="fas fa-globe"></i> Mission Statement</h2>
    <p style='font-size:1.1rem; line-height:1.6; margin-bottom:1.5rem;'>
    Animal testing is still widely used in drug development; but I believe technology can provide better, 
    faster, and more humane alternatives.
    </p>
    <div class="grid-container">
    <div class="tech-card">
    <h4><i class="fas fa-paw"></i> Reduce Animal Testing</h4>
    <p>Minimize need for lab animals in early-stage testing</p>
    </div>
    <div class="tech-card">
    <h4><i class="fas fa-flask"></i> Empower Researchers</h4>
    <p>Provide accessible computational tools for students and scientists</p>
    </div>
    <div class="tech-card">
    <h4><i class="fas fa-rocket"></i> Future of Drug Design</h4>
    <p>Support transition to fully computational approaches</p>
    </div>
    </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Tech Stack Section
    st.markdown("""
    <div class='section'>
    <h2><i class="fas fa-tools"></i> Tech Stack</h2>
    <div class="grid-container">
    <div class="tech-card">
    <h4><i class="fab fa-python"></i> Core</h4>
    <p>Python<br>Scikit-learn<br>RDKit<br>PaDEL</p>
    </div>
    <div class="tech-card">
    <h4><i class="fas fa-desktop"></i> Interface</h4>
    <p>Streamlit<br>HTML/CSS</p>
    </div>
    <div class="tech-card">
    <h4><i class="fas fa-database"></i> Data</h4>
    <p>PubChem<br>Molecular fingerprints<br>SMILES</p>
    </div>
    </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Contact Section
    st.markdown("""
    <div class='section'>
    <h2><i class="fas fa-envelope"></i> Let's Connect!</h2>
    <p style='font-size:1.1rem; text-align:center; margin-bottom:2rem;'>
    I'm always open to collaboration or conversation about bioinformatics, machine learning, or ethical tech.
    </p>
    <div style='display:flex; justify-content:center; gap:3rem; flex-wrap:wrap;'>
    <a href="https://linkedin.com/in/yourprofile" target="_blank" style='text-decoration:none; color:inherit;'>
    <div class="contact-card" style='width:200px;'>
    <i class="fab fa-linkedin social-icon"></i>
    <h4>LinkedIn</h4>
    </div>
    </a>
    <a href="https://github.com/yourprofile" target="_blank" style='text-decoration:none; color:inherit;'>
    <div class="contact-card" style='width:200px;'>
    <i class="fab fa-github social-icon"></i>
    <h4>GitHub</h4>
    </div>
    </a>
    <a href="mailto:your@email.com" style='text-decoration:none; color:inherit;'>
    <div class="contact-card" style='width:200px;'>
    <i class="fas fa-envelope social-icon"></i>
    <h4>Email</h4>
    </div>
    </a>
    </div>
    <p style='text-align:center; margin-top:2rem; font-style:italic; color:#4a8cff;'>
    <i class="fas fa-heart"></i> Thank you for checking out my work! <i class="fas fa-dna"></i>
    </p>
    </div>
    """, unsafe_allow_html=True)

def image_to_base64(image):
    from io import BytesIO
    import base64
    buffered = BytesIO()
    image.save(buffered, format="JPEG")
    return base64.b64encode(buffered.getvalue()).decode()

if __name__ == "__main__":
    about_page()