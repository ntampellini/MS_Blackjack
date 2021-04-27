import streamlit as st
from MS_Blackjack import blackjack_alg

if __name__ == "__main__":

    st.title('MS Blackjack')
    st.write('## Mass Fragments Brute Force Finder')
    st.write('### Nicolò Tampellini - Alessandro Brusa')

    input_string = str(st.number_input('Molecular Peak Mass',value=56))

    st.write('Arguments (optional):')

    if st.checkbox('Hint: uses provided fragment as starting point for formula search. Ex. C6HF2'):
        hint = st.text_area('')
    else:
        hint = False

    exact = st.checkbox('Exact: request exact mass peak (single isotope, best match reported with associate error)')

    mol = st.checkbox(('Mol: requests only structures with integer saturation index (molecules or OE(.+) ions, not EE(+)' + 
                        ' ions) and with a C/H ratio compatible with a proper molecule. Be careful, may discard low-hydrogen' + 
                        'molecules like HC3N or C4N2, which may be proper for human standards.'))

    ext = st.checkbox('Ext: broadens search criterias (required for low %C(m/m) molecules like CS2 and Freons)')

    wild = st.checkbox('Wild: explores all possible chemical space, no assumptions made. (SLOW, USE WITH CARE!)')

    if st.checkbox('Only: use only selected atoms to build structure. Ex. CHClBr'):
        only = st.text_area('')
    else:
        only = False

    nist = st.checkbox('NIST: for each found structure, look up on NIST if there is any matching mass spectra')

    if exact:
        input_string += ' -exact'
    if mol:
        input_string += ' -mol'
    if ext:
        input_string += ' -ext'
    if wild:
        input_string += ' -wild'
    if only:
        input_string += f' -only {only}'
    if hint:
        input_string += f' -hint {hint}'
    if nist:
        input_string += ' -nist'

    if st.button('Run'):
        st.write('Running...')
        output = blackjack_alg(input_string)
        st.write('### Results')
        st.write('\n\n'.join(output))