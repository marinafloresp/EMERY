import seaborn as sns
import streamlit as st
import webbrowser
import altair as alt
import numpy as np
import pandas as pd
import dc_stat_think as dcst

def colours(num):
    if num > 2:
        palette = sns.color_palette("colorblind", num).as_hex()
    else:
        palette = ['#e50000','#02c14d']
    return palette

def add_logo():
    st.markdown(
        """
        <style>
            [data-testid="stSidebarNav"]::before {
                content: "pyBasket";
                margin-left: 20px;
                margin-top: 20px;
                font-size: 30px;
                position: relative;
                top: 50px;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )


def hideRows():
    hide_table_row_index = """
                <style>
                thead tr th:first-child {display:none}
                tbody th {display:none}
                </style>
                """
    return hide_table_row_index

def openGeneCard(gene):
    webpage_link = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + gene
    webbrowser.open(webpage_link)

def savePlot(fig, feature):
    if st.button('Save Plot', key="plot"+feature):
        fig.save('plot_' + feature + '.png')
        st.info('Plot saved as .png in working directory', icon="ℹ️")
    else:
        st.write("")

def savePlot_plt(fig,feature):
    if st.button('Save Plot', key="plot"+feature):
        fig.savefig(feature)
        st.info('Plot saved as .png in working directory', icon="ℹ️")
    else:
        st.write("")

def saveTable(df, feature):
    if st.button('Save table', key="table_"+feature):  # Update the key to a unique value
        df.to_csv('raw_data_'+feature+'.csv', index=False)
        st.info('Data saved as .csv file in working directory', icon="ℹ️")
    else:
        st.write("")

def sideBar():
    if "data" in st.session_state:
        data = st.session_state["data"]
        analysis_data = st.session_state["Analysis"]
        st.sidebar.title("Select basket*cluster interaction")
        with st.sidebar:
            cluster = st.selectbox("Select a cluster", data.clusters_names, key="cluster")
            if "cluster" not in st.session_state:
                st.session_state["cluster"] = cluster
            basket = st.selectbox("Select a basket", data.baskets_names, key="basket")
            if "basket" not in st.session_state:
                st.session_state["basket"] = basket
            subgroup, size = analysis_data.findInteraction(cluster, basket)
            st.info("###### Samples in **cluster {}** & **{} basket**: {}".format(cluster, basket, size))
    else:
        with st.sidebar:
            st.warning("Results not found. Please upload in a results file in Home.")

def alt_line_chart(df, df2, x, y, title_x, title_y, main_title,save):
    base = alt.Chart(df).mark_line().encode(
        alt.X(x+':Q', title=title_x),
        alt.Y(y+':Q', title=title_y)
    ).properties(
        height=650,
        title=main_title,
    )

    percentiles =alt.Chart(df2).mark_point(filled=True, size=300).encode(
        alt.X('x' + ':Q', title=title_x),
        alt.Y('y' + ':Q', title=title_y),color=alt.value('red')
        )
    labels = percentiles.mark_text(
        align='left',
        baseline='middle',
        dx=25,
        fontSize=15
    ).encode(
        text='y'
    )
    base = base + percentiles + labels
    savePlot(base, save)
    st.altair_chart(base, theme="streamlit", use_container_width=True)

def alt_hor_barplot(df, x, y, num_cols, title_x, title_y, colors,main_title,save):
    palette = colours(num_cols)
    base = alt.Chart(df).mark_bar().encode(
        alt.X(x, title=title_x),
        alt.Y(y+':N', title=title_y).sort('-x'), alt.Color(colors +':N')
    ).properties(
        height=650,
        title= main_title
    ).configure_range(
        category=alt.RangeScheme(palette))
    savePlot(base,save)
    st.altair_chart(base, theme="streamlit", use_container_width=True)

def alt_ver_barplot(df, x, y, num_cols, title_x, title_y, colors,main_title,save, tooltip):
    palette = colours(num_cols)
    base = alt.Chart(df).mark_bar().encode(
        alt.X(x +':N', title=title_x),
        alt.Y(y+':Q', title=title_y,axis=alt.Axis(grid=False)), alt.Color(colors+':N'), tooltip = tooltip
    ).properties(
        height=650,
        title=main_title
    ).configure_range(
        category=alt.RangeScheme(palette))
    savePlot(base, save)
    st.altair_chart(base, theme="streamlit", use_container_width=True)

def alt_scatterplot(df, x, y, title_x, title_y,main_title,save,tooltip ):
    base = alt.Chart(df).mark_circle(size=100).encode(
        x=alt.Y(x, title=title_x),
        y=alt.Y(y+':Q', title=title_y),
        tooltip=tooltip
    ).properties(
        height=650,
        title=main_title
    ).interactive().properties(height=650)
    savePlot(base, save)
    st.altair_chart(base, theme="streamlit", use_container_width=True)

def alt_boxplot(df, x, y, num_cols, title_x, title_y, colors,main_title,save):
    palette = colours(num_cols)
    base = alt.Chart(df, title="AAC response").mark_boxplot(extent='min-max', ticks=True, size=60).encode(
        x=alt.X(x, title=title_x),
        y=alt.Y(y, title=title_y), color=alt.Color(colors + ':N')
    ).properties(height=650, title = main_title).configure_range(category=alt.RangeScheme(palette))
    savePlot(base, save)
    st.altair_chart(base, theme="streamlit", use_container_width=True)

def searchTranscripts(transcripts):
    transcript = st.selectbox("Select feature/transcript", transcripts, key="transcript")
    st.caption("Transcripts ordered by decreasing significance.")
    st.write(" ")
    st.write("Click button to search for feature {} in GeneCards database.".format(transcript))
    st.button('Open GeneCards', on_click=openGeneCard, args=(transcript,))
    st.write(" ")
    return transcript

def ecdf(data,intervals):
    percentiles = intervals
    pct_val = np.percentile(data, percentiles)
    x, y = dcst.ecdf(data)
    pct = pd.DataFrame({'Probability': x, 'Percent': y * 100})
    pct_val = pd.DataFrame({'x': np.round(pct_val,3), 'y': percentiles})
    return pct, pct_val
