import streamlit as st
import base64
import os

def get_base64_image(img_path):
    try:
        with open(img_path, "rb") as f:
            return base64.b64encode(f.read()).decode()
    except Exception as e:
        st.warning(f"⚠️ Could not load logo image: {e}")
        return None

def inject_navbar():
    logo_path = "vaccine.png"  # Ensure the path is correct
    logo_data = get_base64_image(logo_path)

    # Remove the default white background and padding
    st.markdown("""
<style>
    /* Adjust Streamlit container width to be full screen */
    .block-container {
        padding-top: 1rem;
        padding-left: 2rem;
        padding-right: 2rem;
        max-width: 100% !important;
    }

    /* Navbar Style */
    .navbar {
        display: flex;
        justify-content: space-between;
        align-items: center;
        background-color: transparent;
        padding: 2rem;
        border: none;
        margin-bottom: -2rem; /* Reduce the gap below the navbar */
    }

    /* Left section: logo and title */
    .navbar-left {
        display: flex;
        align-items: center;
    }

    /* Logo styles */
    .navbar-logo {
        height: 50px;
        margin-right: 5px; /* Reduced margin to bring the logo closer to the text */
        filter: invert(39%) sepia(21%) saturate(469%) hue-rotate(139deg) brightness(92%) contrast(85%);
    }

    /* Title of navbar */
    .navbar-title {
        font-size: 3rem;
        font-weight: 700;
        color: #057a6f;
        margin: 0;
    }

    /* Right section: navbar items */
    .nav-links {
        display: flex;
        gap: 2rem;
    }

    /* Navbar items */
    .nav-item {
        color: #057a6f;
        font-size: 1.2rem;
        font-weight: 600;
        cursor: pointer;
        position: relative;
        transition: color 0.3s ease, transform 0.3s ease;
    }

    /* Hover effect */
    .nav-item:hover {
        color: #045b52;
        text-decoration: none;
        transform: scale(1.1);
    }

    /* Cool underline effect */
    .nav-item::after {
        content: '';
        position: absolute;
        width: 0;
        height: 2px;
        background: linear-gradient(90deg, #057a6f, #045b52);
        bottom: -5px;
        left: 0;
        transition: width 0.4s ease, background-color 0.3s ease;
    }

    .nav-item:hover::after {
        width: 100%;
        background: linear-gradient(90deg, #045b52, #057a6f);
    }
</style>
""", unsafe_allow_html=True)

    # Add logo and title to navbar
    logo_html = f'<img class="navbar-logo" src="data:image/png;base64,{logo_data}">' if logo_data else ''

    # Navbar HTML
    navbar_html = f"""
    <div class="navbar">
        <div class="navbar-left">
            {logo_html}
            <div class="navbar-title">CYPredict</div>
        </div>
        <div class="nav-links">
            <div class="nav-item">Home</div>
            <div class="nav-item">Predictions</div>
            <div class="nav-item">Models</div>
            <div class="nav-item">FAQ</div>
            <div class="nav-item">About Us</div>
            <div class="nav-item">Contact</div>
        </div>
    </div>
    """
    # Inject HTML into Streamlit
    st.markdown(navbar_html, unsafe_allow_html=True)

# Test the navbar


