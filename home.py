import streamlit as st
from streamlit_lottie import st_lottie
import requests
from PIL import Image
import io
import streamlit.components.v1 as components    
import base64
from app import show_navbar


def get_base64_image(img_path):
    with open(img_path, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()


# Set page configuration
st.set_page_config(
    page_title="CYP Predictor",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="collapsed",
)

selected = show_navbar()

# show_navbar() function is not defined. Uncomment and define it if needed.
# show_navbar()
# Custom CSS with card backgrounds using the new color scheme
st.markdown("""
<style>
    @keyframes gradient {
        0% { background-position: 0% 50%; }
        50% { background-position: 100% 50%; }
        100% { background-position: 0% 50%; }
    }
    
    .stat-card {
        border-radius: 12px;
        padding: 1.5rem;
        color: white;
        text-align: center;
        background-size: cover;
        background-position: center;
        position: relative;
        overflow: hidden;
        height: 100%;
        margin-bottom: 2rem; 
    }
            
    .stat-card::before {
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background: rgba(107, 138, 71, 0.7); /* AQUA with 70% opacity */
        z-index: 1;
    }
    
    .stat-content {
        position: relative;
        z-index: 2;
    }
    
    .welfare-card {
        border-radius: 12px;
        padding: 4rem;
        color: white;
        margin: 2rem 0; /* Reduced margin to minimize gap */
        background: linear-gradient(rgba(108, 100, 139, 0.8), rgba(0, 51, 102, 0.8)), /* LAVENDER and DARK BLUE gradient */
                    url('https://aldf.org/wp-content/uploads/2018/05/rabbit-160495596-16x9.jpg');
        background-size: cover;
        background-position: center;
    }
    
    .feature-card {
        border-radius: 12px;
        padding: 1.5rem;
        color: white;
        background-size: cover;
        background-position: center;
        position: relative;
        overflow: hidden;
        height: 100%;
        transition: all 0.3s ease;
    }
    
    .feature-card::before {
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background: rgba(107, 138, 71, 0.7); /* AQUA with 70% opacity */
        z-index: 1;
        transition: all 0.3s ease;
    }
    
    .feature-card:hover::before {
        background: rgba(107, 138, 71, 0.5); /* AQUA with 50% opacity */
    }
    
    .feature-content {
        position: relative;
        z-index: 2;
    }
    
    /* New color scheme classes */
    .aqua-bg {
        background-color: #057a6f;
    }
    
    .lavender-bg {
        background-color: #6C648B;
    }
    
    .darkblue-bg {
        background-color: #003366;
    }
    
    /* Text colors */
    .aqua-text {
        color: #057a6f;
    }
    
    .lavender-text {
        color: #6C648B;
    }
    
    .darkblue-text {
        color: #003366;
    }
</style>
""", unsafe_allow_html=True)

# Hero Section
components.html(open("animated_header.html", "r").read(), height=400)

# Stats Section with Background Images
st.markdown("""
### <span style="font-size: 2rem;">Why Choose CYP Predictor Pro?</span>
""", unsafe_allow_html=True)
cols = st.columns(4)
stats = [
    {
        "value": "80+%", 
        "label": "Prediction Accuracy",
        "bg": "https://img.loigiaihay.com/picture/2024/0729/39-accuracy.jpg",
        "color": "lavender"  # Using LAVENDER for this card
    },
    {
        "value": "30,000+", 
        "label": "Molecular Descriptors Analyzed",
        "bg": "https://images.unsplash.com/photo-1532094349884-543bc11b234d?ixlib=rb-1.2.1&auto=format&fit=crop&w=800&q=80",
        "color": "aqua"  # Using AQUA for this card
    },
    {
        "value": "<2s", 
        "label": "Average Prediction Time per Compound",
        "bg": "https://media.istockphoto.com/id/183876874/photo/hand-with-classic-stopwatch.jpg?s=612x612&w=0&k=20&c=lcebeVhj6sYwHZ6iBM_ph75sUPkb70E1HMv-GhguDMs=",
        "color": "darkblue"  },
    {
        "value": "10+",	 
        "label": "ML Algorithms Tested",
        "bg": "https://images.unsplash.com/photo-1551288049-bebda4e38f71?ixlib=rb-1.2.1&auto=format&fit=crop&w=800&q=80",
        "color": "darkblue"  # Using SUNSHINE for this card
    }
]

for i, stat in enumerate(stats):
    with cols[i]:
        color_map = {
            "aqua": "#057a6f",
            "lavender": "#6C648B",
            "darkblue": "#003366"
        }
        color = color_map[stat["color"]]
        
        st.markdown(f"""
        <div class="stat-card" style="background-image: url('{stat["bg"]}');">
            <div class="stat-content" style="display: flex; flex-direction: column; justify-content: center; height: 100%;">
                <h2 style="margin: 0; font-size: 2.5rem; line-height: 1.2;">{stat['value']}</h2>
                <p style="margin: 1rem 0 0 0; font-size: 1.2rem; line-height: 1.5;">{stat['label']}</p>
            </div>
        </div>
        <style>
            .stat-card:nth-child({i+1})::before {{
                background: rgba({int(color[1:3], 16)}, {int(color[3:5], 16)}, {int(color[5:7], 16)}, 0.7);
            }}
        </style>
        """, unsafe_allow_html=True)

# Animal Welfare Section with Background
try:
    rabbit_base64 = get_base64_image('images/rabbit.png')
    welfare_html = f"""
    <div class="welfare-card">
        <div style="display: flex; align-items: center; justify-content: space-between;">
            <div style="flex: 2; display: flex; flex-direction: column; justify-content: center;">
                <h2 style="color: white; margin-bottom: 2rem; font-size: 2rem; line-height: 1.5;">🐾 Saving Animal Lives Through Innovation</h2>
                <p style="color: white; margin-bottom: 1.5rem; font-size: 1.2rem; line-height: 1.8;">
                    CYPredict is built to reduce reliance on animal testing in drug development by providing accurate, 
                    AI-driven predictions for liver enzyme interactions. Our mission is to empower students, educators, and researchers with tools that prioritize safety, accessibility, and compassion.
                </p>
                <div style="display: flex; gap: 1rem;">
                    <div style="background: rgba(255,255,255,0.2); padding: 1rem; border-radius: 8px; flex: 1;">
                        <p style="margin: 0; font-weight: bold; font-size: 1.1rem;">Reinventing early-stage drug testing</p>
                        <p style="margin: 0; font-size: 1rem; line-height: 1.5;">With data, not animals</p>
                    </div>
                    <div style="background: rgba(255,255,255,0.2); padding: 1rem; border-radius: 8px; flex: 1; ">
                        <p style="margin: 0; font-weight: bold; font-size: 1.1rem;">A learning lab for students, educators, and future pharmacologists</p>
                        <p style="margin: 0; font-size: 1rem; line-height: 1.5;">where science meets empathy</p>
                    </div>
                </div>
            </div>
            <div style="flex: 1; display: flex; justify-content: center;">
                <img src="data:image/png;base64,{rabbit_base64}" width="300">
            </div>
        </div>
    </div>
    """
    st.markdown(welfare_html, unsafe_allow_html=True)
except Exception as e:
    st.error(f"Could not load rabbit image: {str(e)}")
    # Fallback to emoji
    st.markdown("""
    <div class="welfare-card">
        <h2 style="color: white; margin-bottom: 1rem; font-size: 2rem; line-height: 1.5;">🐾 Saving Animal Lives Through Innovation 🐇</h2>
        <!-- rest of the content -->
    </div>
    """, unsafe_allow_html=True)

# Features Section with Background Images
st.markdown("## ✨ Key Features")
features = [
    {
        "title": "Advanced Visualization",
        "description": "Interactive 2D/3D molecule viewer with rotation and zoom capabilities",
        "icon": "🧬",
        "bg": "https://proventainternational.com/wp-content/uploads/2023/08/a.lionsheart_Abstract_3d_render_of_different_molecules_relevant_63450679-4266-43d7-a85c-b64757216f26-1024x683.png",
        "color": "darkblue"  # Using DUSTYROSE for this card
    },
    {
        "title": "Multi-Model Prediction",
        "description": "Predict inhibition for all major CYP450 enzymes simultaneously",
        "icon": "🤖",
        "bg": "https://media.springernature.com/lw900/springer-cms/rest/v1/content/26808036/data/v3",
        "color": "lavender"  # Using LAVENDER for this card
    },
    {
        "title": "Batch Processing",
        "description": "Upload files with multiple compounds for high-throughput screening",
        "icon": "📊",
        "bg": "https://images.unsplash.com/photo-1551288049-bebda4e38f71?ixlib=rb-1.2.1&auto=format&fit=crop&w=800&q=80",
        "color": "darkblue"  # Using SUNSHINE for this card
    }
]

cols = st.columns(3)
for i, feature in enumerate(features):
    with cols[i]:
        st.markdown(f"""
        <div class="feature-card" style="background-image: url('{feature["bg"]}');">
            <div class="feature-content">
                <div style="font-size: 2rem; margin-bottom: 1rem;">{feature['icon']}</div>
                <h3 style="margin-bottom: 0.5rem;">{feature['title']}</h3>
                <p style="line-height: 1.5;">{feature['description']}</p>
            </div>
        </div>
        <style>
            .feature-card:nth-child({i+1})::before {{
                background: rgba(0, 0, 0, 0.5); /* Default shade */
            }}
            .feature-card:nth-child({i+1}):hover::before {{
                background: none; /* Remove shade on hover */
            }}
        </style>
        """, unsafe_allow_html=True)
