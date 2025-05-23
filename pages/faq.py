import streamlit as st
from streamlit_extras.let_it_rain import rain  # type: ignore
from app import show_navbar



# Set page configuration
st.set_page_config(
    page_title="CYP Predictor Pro - FAQ",
    page_icon="❓",
    layout="wide"
)

# Verify correct page is loaded

# Then your navbar and content
selected = show_navbar()

# Custom CSS
st.markdown("""
<style>
    .faq-container {
        max-width: 900px;
        margin: 0 auto;
    }
    
    .faq-question {
        background-color: #057a6f;
        color: white;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
        font-weight: bold;
        cursor: pointer;
        transition: all 0.3s ease;
    }
    
    .faq-question:hover {
        background-color: #04665e;
    }
    
    .faq-answer {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 8px;
        margin-bottom: 1rem;
        border-left: 4px solid #057a6f;
    }
    
    .header-section {
        text-align: center;
        margin-bottom: 2rem;
    }
    
    .emoji-rain {
        font-size: 1.5rem;
    }
</style>
""", unsafe_allow_html=True)

# Header with fun animation
st.markdown("""
<div class="header-section">
    <h1>Frequently Asked Questions ❓</h1>
    <p>Find answers to common questions about CYP Predictor Pro</p>
</div>
""", unsafe_allow_html=True)

# Trigger emoji rain (you'll need streamlit-extras installed)
rain(
    emoji="❓",
    font_size=20,
    falling_speed=5,
    animation_length="infinite",
)

# FAQ Content
faqs = [
    {
        "question": "What is CYP Predictor Pro?",
        "answer": "CYP Predictor Pro is an AI-powered tool that predicts cytochrome P450 enzyme inhibition for drug compounds. It helps researchers and students evaluate potential drug-drug interactions early in the development process."
    },
    {
        "question": "How accurate are the predictions?",
        "answer": "Our models achieve over 80% accuracy across major CYP450 isoforms (1A2, 2C9, 2C19, 2D6, and 3A4). Accuracy may vary depending on the chemical space of your compounds."
    },
    {
        "question": "What input formats does the tool accept?",
        "answer": "You can input compounds as:<br>- SMILES strings with molecular structures"
    },
    {
        "question": "Is there a limit to the number of compounds I can analyze?",
        "answer": "The free version allows up to 10 compounds per batch. Premium accounts can process up to 1,000 compounds at once."
    },
    {
        "question": "What machine learning models are used?",
        "answer": "We employ an ensemble of models including:<br>- Random Forest<br>- XGBoost<br>- SVM <br>Each model is optimized for specific CYP isoforms."
    },
    {
        "question": "How does this tool reduce animal testing?",
        "answer": "By providing accurate in silico predictions, researchers can prioritize compounds for testing, reducing the number of animals needed for preliminary screening by up to 70%."
    },
    {
        "question": "Can I use this for educational purposes?",
        "answer": "Absolutely! We offer special academic licenses and classroom packages. Contact us at education@cyppredictor.com for more information."
    },
    {
        "question": "How do I contact support?",
        "answer": "For any questions or issues, please reach out to our support team. Info below"
    }
]

# Display FAQs
st.markdown('<div class="faq-container">', unsafe_allow_html=True)

for faq in faqs:
    with st.expander(faq["question"]):
        st.markdown(faq["answer"], unsafe_allow_html=True)

st.markdown('</div>', unsafe_allow_html=True)

# Contact section
st.markdown("""
---
### Still have questions?

<div style="background-color:#fcecc0; padding: 1.5rem; border-radius: 8px;">
    <p>Contact our support team:</p>
    <p>📧 <strong>Email:</strong> support@cyppredictor.com</p>
    <p>🐦 <strong>Twitter:</strong> @CYPPredictor</p>
    <p>💻 <strong>Documentation:</strong> <a href="https://docs.cyppredictor.com" target="_blank">docs.cyppredictor.com</a></p>
</div>
""", unsafe_allow_html=True)