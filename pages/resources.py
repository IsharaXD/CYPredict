import streamlit as st
from app import show_navbar


# ========== PAGE CONFIG ==========
st.set_page_config(layout="wide")  # Wide layout



# ========== RESOURCE DATABASE ==========
ARTICLES = [
    {
        "title": "Comprehensive Guide to CYP450 Enzymes",
        "description": "Detailed examination of cytochrome P450 enzyme functions.",
        "tags": ["Metabolism", "Pharmacogenomics"],
        "link": "https://example.com/article1",
        "citation": "Smith et al. (2023) Journal of Pharmacology"
    },
    {
        "title": "Cytochrome P450 Enzymes: Biochemistry, Physiology, and Role in Drug Metabolism",
        "description": "A foundational overview of the P450 superfamily and its pharmacokinetic significance.",
        "tags": ["Biochemistry", "Drug Metabolism"],
        "link": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7371961/",
        "citation": "Zanger & Schwab (2013) Frontiers in Pharmacology"
    },
    {
        "title": "In Vitro Methods for CYP450 Enzyme Induction and Inhibition",
        "description": "Explores experimental models used to study CYP enzyme regulation.",
        "tags": ["In Vitro", "Mechanistic"],
        "link": "https://www.tandfonline.com/doi/abs/10.1080/00498250701534893",
        "citation": "Pelkonen et al. (2008) British Journal of Clinical Pharmacology"
    },
    {
        "title": "Prediction of Human Drug Clearance Using In Silico CYP450 Models",
        "description": "Details machine learning models for predicting drug metabolism routes.",
        "tags": ["In Silico", "Pharmacokinetics", "Modeling"],
        "link": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5696450/",
        "citation": "Kirchmair et al. (2015) Computational and Structural Biotechnology Journal"
    }
]

VIDEOS = [
    {
        "title": "CYP450 Enzymes Drug Interactions MADE EASY in 5 MINS",
        "description": "Animated overview of CYP enzyme functions (15 min).",
        "tags": ["Education", "Animation"],
        "youtube_id": "vle_0dN3bwA",
        "duration": "5:16"
    },
    {
        "title": "Pharmacology - DRUG INTERACTIONS (MADE EASY)",
        "description": "Simple explanation of CYP isoenzymes, families, and clinical importance.",
        "tags": ["Education", "Pharmacology"],
        "youtube_id": "BayzMj3n5_I",
        "duration": "13:30"
    },
    
    {
        "title": "In Silico Metabolism Webinar",
        "description": "The Metabolism Module in ADMET Predictor™ contains in silico models that classify compounds as substrates and/or inhibitors of the major CYP isoforms, while also predicting likely sites of metabolism and kinetic parameters (Km, Vmax, and intrinsic clearance). Classification models for phase II glucuronidation by UDP-glucuronosyltransferase are also included. This webinar describes the development of these models and how they can be applied to predict the disposition of drug candidates and assist with the lead optimization process.",
        "tags": ["AI", "Modeling", "CYP Inhibition"],
        "youtube_id": "BT9CzCsi6lU",
        "duration": "1:04:15"
    }
]

# ========== STYLES ==========
st.markdown("""
<style>
    /* Main container */
    .main-container {
        padding: 2rem 5%;
        max-width: 1400px;
        margin: 0 auto;
    }
    
    /* Header animation */
    @keyframes fadeInSlide {
        0% { opacity: 0; transform: translateY(-30px); }
        50% { opacity: 0.5; transform: translateY(-15px); }
        100% { opacity: 1; transform: translateY(0); }
    }
    .animated-header {
        animation: fadeInSlide 1.5s ease-out;
        text-align: center; /* Center the header */
    }
    
    
    /* Resource cards */
    .resource-card {
        border-radius: 12px;
        padding: 2rem;
        margin: 1.5rem 0;
        background: white;
        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        transition: transform 0.4s ease, box-shadow 0.4s ease;
    }
    .resource-card:hover {
        transform: translateY(-5px) scale(1.02);
        box-shadow: 0 6px 16px rgba(0,0,0,0.15);
    }
    
    /* Titles */
    .card-title {
        font-size: 1.5rem;
        color: #057a6f;
        margin-bottom: 1rem;
        font-weight: 700;
    }
    
    /* Descriptions */
    .card-description {
        font-size: 1.1rem;
        line-height: 1.6;
        color: #333;
        margin-bottom: 1.2rem;
    }
    
    /* Tags */
    .tag {
        display: inline-block;
        background: #e6f7f5;
        color: #057a6f;
        padding: 0.4rem 1rem;
        border-radius: 20px;
        margin-right: 0.8rem;
        font-size: 0.9rem;
    }
    
    /* Action buttons */
    .action-btn {
        display: inline-block;
        background: #057a6f;
        color: white !important;
        padding: 0.6rem 1.5rem;
        border-radius: 8px;
        text-decoration: none !important;
        margin-top: 1rem;
        font-weight: 600;
        transition: background 0.3s, transform 0.3s;
    }
    .action-btn:hover {
        background: #04665e;
        transform: scale(1.05);
    }
    
    /* Video container */
    .video-container {
        margin: 1rem 0;
        border-radius: 10px;
        overflow: hidden;
        position: relative;
        padding-bottom: 40%;
    }
    .video-container iframe {
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
    }
    
</style>
""", unsafe_allow_html=True)

# ========== HELPER FUNCTIONS ==========
def display_article(article):
    """Show an article card"""
    st.markdown(f"""
    <div class="resource-card">
        <div class="card-title">{article['title']}</div>
        <div class="card-description">{article['description']}</div>
        <div>
            {"".join(f"<span class='tag'>{tag}</span>" for tag in article['tags'])}
        </div>
        {f"<div style='margin-top: 1rem; color: #666;'><i>{article['citation']}</i></div>" 
         if article.get('citation') else ""}
        <a href="{article['link']}" class="action-btn" target="_blank">
           🌐 Read Article
        </a>
    </div>
    """, unsafe_allow_html=True)

def display_video(video):
    """Show a video card"""
    st.markdown(f"""
    <div class="resource-card">
        <div class="card-title">{video['title']}</div>
        <div class="card-description">
            {video['description']}
            {f"<div style='margin-top: 0.5rem; color: #666;'>Duration: {video['duration']}</div>" 
             if video.get('duration') else ""}
        </div>
        <div>
            {"".join(f"<span class='tag'>{tag}</span>" for tag in video['tags'])}
        </div>
        <div class="video-container">
            <iframe src="https://www.youtube.com/embed/{video['youtube_id']}" 
                    frameborder="0" allowfullscreen>
            </iframe>
        </div>
    </div>
    """, unsafe_allow_html=True)

# ========== MAIN PAGE ==========
def show_resources():
    show_navbar()
    
    st.markdown("""
    <div class="main-container">
        <div class="animated-header">
            <h1 style='color: #057a6f; margin-bottom: 0.5rem;'>📚 Educational Resources</h1>
            <p style='font-size: 1.1rem; color: #555; margin-bottom: 2rem;'>
                Curated collection of learning materials about CYP enzymes
            </p>
        </div>
    """, unsafe_allow_html=True)
    
    # Search and filters
    search_query = st.text_input("🔍 Search resources", placeholder="Search by title, tags...")
    selected_tags = st.multiselect(
        "🏷 Filter by tags",
        options=sorted(list(set(tag for item in ARTICLES + VIDEOS for tag in item["tags"])))
    )
    
    # Tabs
    tab1, tab2 = st.tabs(["📖 Articles", "🎥 Videos"])
    
    with tab1:
        filtered_articles = [
            a for a in ARTICLES 
            if (not search_query or search_query.lower() in a["title"].lower()) and
               (not selected_tags or any(tag in selected_tags for tag in a["tags"]))
        ]
        
        if not filtered_articles:
            st.info("No articles found matching your criteria")
        for article in filtered_articles:
            display_article(article)
    
    with tab2:
        filtered_videos = [
            v for v in VIDEOS 
            if (not search_query or search_query.lower() in v["title"].lower()) and
               (not selected_tags or any(tag in selected_tags for tag in v["tags"]))
        ]
        
        if not filtered_videos:
            st.info("No videos found matching your criteria")
        for video in filtered_videos:
            display_video(video)
    
    st.markdown("</div>", unsafe_allow_html=True)

if __name__ == "__main__":
    show_resources()