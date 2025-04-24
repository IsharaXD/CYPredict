import streamlit as st
import base64

def get_base64_image(img_path):
    try:
        with open(img_path, "rb") as f:
            return base64.b64encode(f.read()).decode()
    except Exception as e:
        return None

def show_navbar():
    logo_data = get_base64_image('images/vaccine.png')
    
    st.markdown(f"""
    <style>
        /* Remove all default Streamlit UI */
        #MainMenu, header, footer {{ visibility: hidden; }}
        .stApp [data-testid="stToolbar"] {{ display: none; }}
        
        /* Completely remove sidebar and its toggle */
        section[data-testid="stSidebar"], .stApp [data-testid="collapsedControl"] {{
            display: none !important;
        }}
        
        /* Reset app padding */
        .stApp {{
            padding-top: 0 !important;
            margin-top: 0 !important;
        }}
        
        /* Navbar container - now scrolls with page */
        .navbar {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-top: -20px;
            padding: 0.8rem 3rem;
            background-color: white;
            position: sticky;  /* Changed from fixed to sticky */
            top: 0;  /* Stick to top when scrolling */
            z-index: 1000;
            border: none !important;
            box-shadow: none !important;
            outline: none !important;
            width: 100%;
        }}
        
        /* Adjust content spacing */
        .stApp > div:first-child {{
            padding-top: 0 !important;
        }}
        
        /* Logo styling */
        .navbar-logo {{
            height: 60px;
            margin-right: -8px;
            filter: invert(39%) sepia(21%) saturate(469%) hue-rotate(139deg) brightness(92%) contrast(85%);
        }}
        
        /* Navigation links */
        .nav-links {{
            display: flex;
            gap: 2rem;
        }}
        
        /* Link styling */
        .nav-link, .nav-link:link, .nav-link:visited, .nav-link:hover, .nav-link:active {{
            color: #057a6f !important;
            font-size: 1.3rem;
            font-weight: 600;
            text-decoration: none !important;
            cursor: pointer;
            padding: 0.5rem 0;
            position: relative;
            border: none !important;
            outline: none !important;
        }}
        
        /* Underline animation */
        .nav-link::after {{
            content: '';
            position: absolute;
            width: 0;
            height: 2px;
            background: #057a6f;
            bottom: 0;
            left: 0;
            transition: width 0.3s;
        }}
        .nav-link:hover::after {{
            width: 100%;
        }}
    </style>

    <!-- Navbar HTML -->
    <div class="navbar">
        <a href="/" style="display: flex; align-items: center; gap: 10px; text-decoration: none; border: none !important;">
            <img class="navbar-logo" src="data:image/png;base64,{logo_data or ''}" style="border: none !important;">
            <span style="font-size: 3rem; font-weight: 700; color: #057a6f; border: none !important;">CYPredict</span>
        </a>
        <div class="nav-links">
            <a href="/" class="nav-link" target="_self">Home</a>
            <a href="/predictions" class="nav-link" target="_self">Predict</a>
            <a href="/Models" class="nav-link" target="_self">Models</a>
            <a href="/resources" class="nav-link" target="_self">Resources</a>
            <a href="/faq" class="nav-link" target="_self">FAQ</a>
        </div>
    </div>
    """, unsafe_allow_html=True)