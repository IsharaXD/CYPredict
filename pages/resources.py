import streamlit as st
from app import show_navbar # Your existing navbar component

# Initialize session state for search/filters
if 'search_query' not in st.session_state:
    st.session_state.search_query = ""
if 'selected_tags' not in st.session_state:
    st.session_state.selected_tags = []

def show_resources():
    # Show your navbar
    show_navbar()
    
    st.markdown("""
    <style>
        /* Custom card styling */
        .resource-card {
            border-radius: 10px;
            padding: 1.5rem;
            margin-bottom: 1rem;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            transition: transform 0.2s;
        }
        .resource-card:hover {
            transform: translateY(-3px);
        }
        .tag {
            display: inline-block;
            background-color: #057a6f;
            color: white;
            padding: 0.2rem 0.6rem;
            border-radius: 20px;
            font-size: 0.8rem;
            margin-right: 0.5rem;
            margin-bottom: 0.5rem;
        }
        .featured-badge {
            background-color: #ffd700;
            color: #333;
        }
    </style>
    """, unsafe_allow_html=True)

    # Page header
    st.title("📚 Resources Center")
    st.markdown("Explore educational materials about CYP enzymes and drug metabolism")

    # Tabs
    tab1, tab2 = st.tabs(["📖 Articles", "🎥 Videos"])

    # Sample data - replace with your actual data
    articles = [
        {
            "title": "Understanding CYP450 Enzymes",
            "description": "Comprehensive guide to cytochrome P450 enzyme functions",
            "tags": ["Metabolism", "Toxicology"],
            "link": "https://example.com/article1",
            "featured": True,
            "type": "pdf"
        },
        # Add more articles...
    ]

    videos = [
        {
            "title": "CYP450 Drug Interactions Explained",
            "description": "Visual guide to common drug interactions",
            "tags": ["Education", "Pharmacology"],
            "youtube_id": "dQw4w9WgXcQ",  # Replace with actual IDs
            "featured": True
        },
        # Add more videos...
    ]

    # Search and filters (applies to both tabs)
    col1, col2 = st.columns([3, 1])
    with col1:
        st.session_state.search_query = st.text_input("🔍 Search resources", placeholder="Search by title, description...")
    with col2:
        all_tags = list(set(tag for item in articles + videos for tag in item["tags"]))
        st.session_state.selected_tags = st.multiselect("🏷 Filter by tags", all_tags)

    with tab1:
        st.header("Educational Articles")
        
        # Filter articles
        filtered_articles = [
            a for a in articles 
            if (st.session_state.search_query.lower() in a["title"].lower() or 
                st.session_state.search_query.lower() in a["description"].lower()) and
               (not st.session_state.selected_tags or 
                any(tag in st.session_state.selected_tags for tag in a["tags"]))
        ]

        for article in filtered_articles:
            with st.container():
                st.markdown(f"""
                <div class="resource-card">
                    <h3>{article['title']}</h3>
                    <p>{article['description']}</p>
                    <div>
                        {"<span class='tag featured-badge'>Featured</span>" if article.get('featured') else ""}
                        {"".join(f"<span class='tag'>{tag}</span>" for tag in article['tags'])}
                    </div>
                    <div style='margin-top: 1rem;'>
                        {"📄 <a href='{article['link']}' target='_blank'>Read Online</a>" if article['type'] == 'link' else ""}
                        {"📥 <a href='{article['link']}' download>Download PDF</a>" if article['type'] == 'pdf' else ""}
                    </div>
                </div>
                """, unsafe_allow_html=True)

    with tab2:
        st.header("Educational Videos")
        
        # Filter videos
        filtered_videos = [
            v for v in videos 
            if (st.session_state.search_query.lower() in v["title"].lower() or 
                st.session_state.search_query.lower() in v["description"].lower()) and
               (not st.session_state.selected_tags or 
                any(tag in st.session_state.selected_tags for tag in v["tags"]))
        ]

        for video in filtered_videos:
            with st.container():
                st.markdown(f"""
                <div class="resource-card">
                    <h3>{video['title']}</h3>
                    <p>{video['description']}</p>
                    <div>
                        {"<span class='tag featured-badge'>Featured</span>" if video.get('featured') else ""}
                        {"".join(f"<span class='tag'>{tag}</span>" for tag in video['tags'])}
                    </div>
                    <div style='margin-top: 1rem;'>
                        <iframe width="100%" height="315" src="https://www.youtube.com/embed/{video['youtube_id']}" 
                        frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" 
                        allowfullscreen></iframe>
                    </div>
                </div>
                """, unsafe_allow_html=True)

if __name__ == "__main__":
    show_resources()