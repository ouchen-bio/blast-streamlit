import streamlit as st
import io
import base64
import pandas as pd
from streamlit_option_menu import option_menu
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from http.client import IncompleteRead

def main():
    
    #Page config
    st.set_page_config(layout='wide')
    
    logo=open("logo.png","rb")
    l=logo.read()
    dat=base64.b64encode(l).decode("UTF-8")
    logo.close()
    st.sidebar.markdown(f'<img src="data:image/png;base64,{dat}" width="100%" alt="logo png">',unsafe_allow_html=True)
    st.sidebar.markdown('')
    st.sidebar.markdown('')
    
    #Menu bar
    with st.sidebar:
        main_menu = option_menu(menu_title='',options=['About','Nucleotide BLAST','Protein BLAST'],
                                icons=['','',''],
                                default_index=0,
                                styles={'container':{'background-color':'#FFFFFF'},
                                                        'nav-link':{'font-size':'14px','text-align':'left'},
                                                        'nav-link-selected':{'background-color':'#0557FD','text_color':'2d2d2d'}})
    
    def blast_search(algorithm, database ,query ,threshold_identity, threshold_coverage):
        try:
            #Perform the BLAST search
            result_handle = NCBIWWW.qblast(algorithm, database, query)
            #Parse the BLAST results
            blast_record = NCBIXML.read(result_handle)
            
        ####### Handle errors starts ########
        except IncompleteRead:
            st.write('Check Your Connexion And Try Again')
        except Exception as e:
            st.write(f"An error occurred: {str(e)}")
        
        ####### Handle errors ends ########
        else:
            filtered_hits = []
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    query_length = len(query)
                    identity=(float(hsp.identities) / float(hsp.align_length)) * 100.0
                    coverage=(float(hsp.align_length) / float(query_length)) * 100.0
                    if identity >= threshold_identity and coverage >= threshold_coverage:
                        filtered_hits.append((alignment.title, hsp.score, hsp.bits, identity, coverage, hsp.expect,
                                                 hsp.align_length, query_length,hsp.query_start,hsp.query_end, hsp.query[0:75] + '...',
                                                 hsp.match[0:75] + '...', hsp.sbjct[0:75] + '...', alignment.length, hsp.sbjct_start, 
                                              hsp.sbjct_end))
            if len(filtered_hits) > 0:
                st.subheader("Search results")
                st.markdown("***")
                for hit in filtered_hits:
                    st.write("**Alignment**: ", hit[0])
                    data = [[hit[1],hit[2],hit[3],hit[4],hit[5],hit[6],hit[7],hit[13]]]
                    statistics = pd.DataFrame(data,columns=["Score","Bit Score","Identity (%)","Coverage (%)","E-Value",
                                                            "Alignment Length (bp)","Query Length (bp)","Subject Length (bp)"])
                    st.dataframe(statistics)
                    st.write("**Query start**: ", hit[8])
                    st.write("**Query end**: ", hit[9])
                    st.write(hit[10])
                    st.write(hit[12])
                    st.write("**Subject start**: ", hit[14])
                    st.write("**Subject end**: ", hit[15])
                    st.markdown("***")
            else:
                st.write("No hits found.")
    
    if main_menu == 'About':
        st.sidebar.info('''
        This app is made using :
        
        Streamlit - Biopython - Pandas
        
        Email: ouchen.yassine.umi@gmail.com
        LinkedIn: [Click Here](https://www.linkedin.com/in/yassine-ouchen-be)
        ''')
        st.header('About BLAST-Streamlit')
        st.write('The BLAST-Streamlit app is a tool for performing sequence similarity searches using the Basic Local Alignment Search Tool (BLAST) algorithm from Biopython. It allows you to search for similarities between nucleotide or protein sequences of interest and sequences in various databases.')
        st.write('For more information about BLAST, please refer to the [NCBI BLAST website](https://blast.ncbi.nlm.nih.gov/Blast.cgi).')
        st.markdown('')
        st.header('Features')
        st.write('- **Nucleotide BLAST**: Perform nucleotide sequence similarity searches using the `blastn` algorithm.')
        st.write('- **Protein BLAST**: Perform protein sequence similarity searches using the `blastp` algorithm.')
        st.write('- **Input Sequence**: Enter a sequence directly into the app for the BLAST search.')
        st.write('- **Import FASTA**: Upload a FASTA file containing a sequence for the BLAST searches.')
        st.write('- **Database Selection**: Choose from a variety of databases to search against, including `nt`, `est`, `TSA` for nucleotide BLAST, and `nr`, `env_nr`, `tsa_nr` for protein BLAST.')
        st.write('- **Thresholds**: Set identity percentage and coverage percentage thresholds to filter the search results.')
        st.write('- **Search Results**: View the alignment details, statistics, and hit sequences for the matching results.')
        st.markdown('')
        st.header('How to Use BLAST-Streamlit')
        st.write('1. Select the type of BLAST search (nucleotide or protein) from the main menu.')
        st.write('2. Choose whether to input the sequence directly or upload a FASTA file.')
        st.write('3. Enter the sequence or upload the file accordingly.')
        st.write('4. Select the database to search against.')
        st.write('5. Adjust the identity and coverage thresholds if desired.')
        st.write('6. Click on the "Search" button to perform the BLAST search.')
        st.write('7. View the search results, including alignment details and statistics.')
     
        
    elif main_menu == 'Nucleotide BLAST':
        algorithm = 'blastn'
        gene_entry = st.sidebar.radio('',["Input Sequence","Import FASTA"],horizontal=True)
        
        if gene_entry=="Input Sequence":
            input_seq = st.sidebar.text_area('',height=120)    
            query = input_seq.upper()
            
            database = st.sidebar.selectbox('Database',['nt','est','TSA'])
            threshold_identity = st.sidebar.slider('Identity %', 50, 100, 80)
            threshold_coverage = st.sidebar.slider('Coverage %', 50, 100, 80)
          
            
            
            if st.sidebar.button('🔍 Search'):
                #Check if the user has entered a sequence
                if not input_seq:
                    st.warning("warning: Please Enter a Gene Sequence.")
                else:
                    blast_search(algorithm, database,query, threshold_identity, threshold_coverage)
                    
            st.sidebar.markdown('')
        else:
            # Create Uploader
            file_uploader=st.sidebar.file_uploader('',type=["fasta"])
            
            database = st.sidebar.selectbox('Database',['nt','est','TSA'])
            threshold_identity = st.sidebar.slider('Identity %', 50, 100, 80)
            threshold_coverage = st.sidebar.slider('Coverage %', 50, 100, 80)
            
            if st.sidebar.button('🔍 Search'):
                if not file_uploader:
                    st.warning("warning: Please Upload a FASTA File.")
                else:
                    #Converting Genbank Format Into A Readable Streamlit Format
                    if file_uploader is not None:
                        byte_str=file_uploader.read()
                        text_obj=byte_str.decode("UTF-8")
                        seq_object=SeqIO.parse(io.StringIO(text_obj),"fasta")
                        sequences=[]
                        for record in seq_object:
                            sequences.append(record.seq)
                        for query in sequences:
                            blast_search(algorithm, database, query,threshold_identity, threshold_coverage)
            st.sidebar.markdown('')
    else:
        algorithm = 'blastp'
        gene_entry = st.sidebar.radio('',["Input Sequence","Import FASTA"],horizontal=True)
        if gene_entry=="Input Sequence":
            input_seq = st.sidebar.text_area('',height=120)    
            query = input_seq.upper()
            
            database = st.sidebar.selectbox('Database',['nr','env_nr','tsa_nr'])
            threshold_identity = st.sidebar.slider('Identity %', 50, 100, 80)
            threshold_coverage = st.sidebar.slider('Coverage %', 50, 100, 80)
            
            
            
            if st.sidebar.button('🔍 Search'):
                #Check if the user has entered a sequence
                if not input_seq:
                    st.warning("warning: Please Enter a Protein Sequence.")
                else:
                    blast_search(algorithm, database, query,threshold_identity, threshold_coverage)
            st.sidebar.markdown('')        
        else:
            # Create Uploader
            file_uploader=st.sidebar.file_uploader('',type=["fasta"])
            
            database = st.sidebar.selectbox('Database',['nr','env_nr','tsa_nr'])
            threshold_identity = st.sidebar.slider('Identity %', 50, 100, 80)
            threshold_coverage = st.sidebar.slider('Coverage %', 50, 100, 80)
            
            if st.sidebar.button('🔍 Search'):
                if not file_uploader:
                    st.warning("warning: Please Upload a FASTA File.")
                else:
                    #Converting Genbank Format Into A Readable Streamlit Format
                    if file_uploader is not None:
                        byte_str=file_uploader.read()
                        text_obj=byte_str.decode("UTF-8")
                        seq_object=SeqIO.parse(io.StringIO(text_obj),"fasta")
                        sequences=[]
                        for record in seq_object:
                            sequences.append(record.seq)
                        for query in sequences:
                            blast_search(algorithm, database, query, threshold_identity, threshold_coverage)
            st.sidebar.markdown('')                
    
if __name__ == '__main__':
    main()
