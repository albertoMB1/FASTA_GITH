
import os
import pandas as pd
import streamlit as st
import pandas as pd
from io import StringIO
from LIBRERIA import fasta1 as fst
import orjson
import re
from LIBRERIA import vcf41 as vc
import time
import orjson
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import altair as alt
import math

# Configuración página inicial
st.set_page_config(
                    page_title="TFG UPV",
                    initial_sidebar_state="expanded",
                    layout='wide'
                    )
# Paginas de la app
PAGES = [
    'Descarga Secuencia Cromosoma',
    'Buscar secuencia',
    'Substituir',
    'K-mer'
]

# Diccionario radio button
dic ={'1': 'Descarga Secuencia Cromosoma',
          '2': 'Buscar secuencia',
          '3': 'Substituir',
          '4':'K-mer'}

def validar_email(email):
    patron = r'^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*\.[a-zA-Z]+$'
    return re.match(patron, email) is not None

@st.cache_data
def leer_orjson(ruta):
    with open(ruta, 'rb') as f:
        dic_indice = orjson.loads(f.read())
    print(f"Diccionario leido en {ruta}")
    return dic_indice


RUTA = 'VCF/GRCh38_latest_clinvar.vcf'
ruta_rango_crom = 'DATOS/Nom_crom/dic_crom_rango.json'
ruta_crom_offline = 'DATOS/GENOMA'

ruta_archivo_clndn = "DATOS/CLNDN/CLNDN.json"
if not os.path.exists(RUTA):
    st.error("Copia el archivo VCF en la carpeta VCF con el nombre GRCh38_latest_clinvar.vcf")

if os.path.exists(ruta_archivo_clndn):
        indices_enfermedades_clndn = leer_orjson(ruta_archivo_clndn)

else:
    st.error("ATENCIÓN: Ejecuta PRIMERO análisis -> CLNDN de app VCF y Copia el diccionario en DATOS\CLNDN")

diccionario_correspondencia = {
        'NC_000007': '7', 'NC_000023': 'X', 'NT_187661.1': 'NT_187661.1', 'NC_000021': '21', 'NC_000005': '5',
        'NC_000001': '1', 'NC_000018': '18', 'NC_000003': '3', 'NC_000020': '20', 'NC_000004': '4','NC_000006': '6',
        'NC_000022': '22', 'NW_009646201.1': 'NW_009646201.1', 'NC_000002': '2','NC_000024': 'Y',
        'NC_000019': '19', 'NC_000015': '15', 'NC_000017': '17', 'NT_187633.1': 'NT_187633.1',
        'NC_012920': 'MT', 'NC_000013': '13','NT_113889.1': 'NT_113889.1','NC_000008': '8',
        'NT_187693.1': 'NT_187693.1', 'NC_000011': '11', 'NC_000016': '16', 'NC_000014': '14',
        'NC_000009': '9', 'NC_000010': '10', 'NC_000012': '12'
    }
columnas_busqueda = ["#CHROM", "#POS", "#ID", "#REF", "#ALT", "#QUAL", "#FILTER", "#INFO"]
# Comprueba si el archivo existe
if os.path.exists(ruta_rango_crom):
    dic_rango_indices = leer_orjson(ruta_rango_crom)
else:
    # El archivo no existe, muestra un mensaje de error
    st.error("Atención: El archivo 'dic_crom_rango.json' no existe: Ejecuta escarga Secuencia Cromosoma -> Obtener Rangos Cromosomas -> De todos los cromosomas")



nombre_cromosomas_vcf = {'NT_187633.1', '22', 'MT', '11', '21', '10', '13', '5', '20', 'X', '3', 'NT_187661.1', 'Y', 'NT_113889.1', '9', '1', 'NT_187693.1', '18', 'NW_009646201.1', '15', '2', '14', '7', '6', '19', '17', '4', '12', '8', '16'}


if os.path.exists(ruta_rango_crom):
    dic_rango_crom = leer_orjson(ruta_rango_crom)
else:
    st.error("ATENCIÓN: Descarga los rangos de cromosomas primero")



def run():

    st.sidebar.title('Luis Alberto MB')
    if st.session_state.page:
        page=st.sidebar.radio('## Navegación', PAGES, index=st.session_state.page)
    else:
        page=st.sidebar.radio('## Navegación', PAGES, index=0,)
    

    if page == 'Descarga Secuencia Cromosoma':
        st.title("Consulta y descarga partes del cromosoma")
        st.sidebar.write("""
            ## Sobre:
            Acciones que se ejecutan:\n
                    """)
        st.sidebar.write("""- Creación Diccionario Rangos Cromosoma. Consulta de Rango Máximo por cromosoma
                        \n- Descarga todos los cromosomas.
                         \n- Crea archivo de cromosoma definidas en un rango.
                         """)
        subcolumna1, subcolumna2 = st.columns(2)
        with subcolumna1:
            st.session_state.data_type = st.radio("Elige una opción", ('Obtener Rangos Cromosomas','Online', 'OffLine'), index = 1)

        if st.session_state.data_type =='Obtener Rangos Cromosomas':
            with subcolumna2:
                st.session_state.data_format = st.radio('Elige Opción', ['De todos los cromosomas','Selecciona Cromosomas'], index = 0)
            if st.session_state.data_format == 'De todos los cromosomas':
                    st.write('---')
                    st.write(f'El email seleccionado: {st.session_state.email}')
                    if st.button('Descargar rangos de cromosomas'):
                        #conj_crom = vc.obtain_all_names_chrom_inVCF(RUTA)
                        dic_rango_aux = crea_muestra_dic_rangos(nombre_cromosomas_vcf)
                        vc.guardar_dicc_json(ruta_rango_crom, dic_rango_aux)
                        st.success(f'Archivo de rango guardado en {ruta_rango_crom}')

            if st.session_state.data_format == 'Selecciona Cromosomas':
                    st.write('---')
                    st.write('##### Selecciona los cromosomas individualmente de los que quieres conocer su rango')
                    seleccionados = st.multiselect('Selección de cromosomas',sorted(nombre_cromosomas_vcf),placeholder='Selecciona Cromosoma')
       
                    if not seleccionados == []:
                        if st.button('Mostrar Seleccionado'):
                            dic_rango_aux = crea_muestra_dic_rangos(seleccionados)
                            diccionario_bytes = orjson.dumps(dic_rango_aux)

                            # Crear un botón de descarga para el archivo json serializado
                            st.download_button(
                                label="Descargar datos como JSON",
                                data=diccionario_bytes,
                                file_name="datos.json",
                                mime="application/json")

        if st.session_state.data_type =='Online':
            with subcolumna2:
                st.session_state.data_format = st.radio('Descargar', ['Todo el Genoma Completo', 'Secuencia Cromosoma'],0)
            if st.session_state.data_format == 'Todo el Genoma Completo':
                    st.write('---')
                    st.toast(f'Esto puede tardar más de 30 minutos', icon='⌛')
                    if st.button('Descargar Genoma'):
                        for clave, valor in dic_rango_indices.items():
                            st.toast(f'Iniciando descarga de cromosoma {clave}', icon='🕠')
                            h = time.time()
                            fst.fetch_genomic_dna_archivo(clave,1, valor,st.session_state.email,f'DATOS/GENOMA/{clave}')
                            j = time.time()
                            st.toast(f'Cromosoma {clave}, descargado en DATOS/GENOMA/{clave}.fasta en {round(j-h,2)} segundos', icon='✅' )
                        st.success(f'Descargado Genoma Completo')


            if st.session_state.email:
                    st.write("Email guardado:", st.session_state.email)
            else:
                    st.write("No se ha proporcionado ningún email.")

            if st.session_state.data_format =='Secuencia Cromosoma':

                st.write('---')
                nom_ncbi = fst.convertir_nombre_crom_del_vcf_a_ncbi(nombre_cromosomas_vcf)
                st.write('#### Elige los cromosomas')
                seleccionados = st.multiselect('Elige cromosomas',sorted(nom_ncbi),placeholder='Selecciona Cromosoma')

                if not seleccionados == []:
                    st.write('''##### Indica el inicio y final del cromosoma que quieres descargar ''')
                    dic_query = seleccion_rangos(seleccionados)

                    if st.button('Consultar Secuencia Cromosomas Seleccionados') and comprobar_seleccion_rangos(dic_query):

                        for clave, valor in dic_query.items():
                            st.toast(f'Consultando rango de cromosoma {clave}', icon='🕠')
                            h = time.time()
                            fst.fetch_genomic_dna_archivo(clave,valor[0], valor[1],st.session_state.email,f'Output/{clave}_{valor[0]}_{valor[1]}')
                            j = time.time()
                            st.success(f'Cromosoma {clave} descargado en {round(j-h,2)} segundos. Desde posición {valor[0]} hasta {valor[1]}')
                            st.toast(f'Descargado Cromosoma {clave} desde la posición {valor[0]}, hasta la posición {valor[1]}',icon='✅' ) 
                            st.success(f'El archivo {clave}_{valor[0]}_{valor[1]}.fasta se ha descargado en la carpeta Output')
                        


        if st.session_state.data_type =='OffLine':
  
            archivos_fasta = listar_archivos_fasta(ruta_crom_offline)
        
            if not archivos_fasta: st.error("No existen archivos fasta descargados previamente. Pulsa online-> Descargar todo el genoma")
            nombres_sin_extension = [clave for clave, _ in archivos_fasta.items()]

            #print(nombres_sin_extension)
            st.write('---')
            st.write('### Elige los cromosomas')
            seleccionados = st.multiselect('Elige Cromosomas', nombres_sin_extension, placeholder='Seleccionar Cromosa')
            st.write('''##### Indica el inicio y final de la secuencia del cromosoma que quieres extraer ''')
            dic_query = seleccion_rangos(seleccionados)
            
            if not seleccionados ==[]:
                if st.button('Consultar subsecuencias'):
                    columnas_seq_extra = ['Cromosoma', 'Inicio', 'Fin', 'Secuencia extraída']
                    df = pd.DataFrame(columns=columnas_seq_extra)
                    
                    st.write('---')
                    for crom_selec, secuencia in dic_query.items():
                        sequens = fst.extraer_subsecuencia_de_fasta_file(archivos_fasta[crom_selec][1],secuencia[0], secuencia[1])
                        df = anyadir_fila(df, crom_selec, secuencia[0], secuencia[1],sequens)
                        guardar_archivo_fasta(crom_selec, secuencia[0], secuencia[1], sequens)
                        st.toast(f'Archivo Ouput/{crom_selec}_{secuencia[0]}_{secuencia[1]}.fasta guardado')
                    st.success(f'Subsecuencias obtenidas')
                    st.dataframe(df)
                    st.success(f'Archivos guardados en carpeta Ouput')
                    








#--------------------------------------BUSCAR SECUENCIA---------------------------------------------
    if page == 'Buscar secuencia':
        st.sidebar.write("""
            ## Sobre:
            Búsqueda de diferentes campos:\n
                    """)
        st.sidebar.write("""
                         \n- Búsqueda de Secuencia Específica en cromosomas
                        \n- De forma En Línea o Local.
                        
                         """)
        st.title("Búsquedas")

        subcolom1, subcolom2 = st.columns(2 )
        with subcolom1:
            st.session_state.data_type = st.radio("Elige una ítem para buscar por:", ('OnLine',  'OffLine'), index = 1)

  
        if st.session_state.data_type =='OnLine':
            st.write('---')
            nom_ncbi = fst.convertir_nombre_crom_del_vcf_a_ncbi(nombre_cromosomas_vcf)
            st.write('#### Elige los cromosomas')

            seleccionados = st.multiselect('Elige cromosomas',sorted(nom_ncbi),placeholder='Selecciona Cromosoma')
            if not seleccionados == []:
                st.write('''##### Indica el inicio y final del cromosoma que quieres consultar ''')
                dic_query = seleccion_rangos(seleccionados)
             
                cadena = st.text_input('Introduce cadena a buscar').upper().replace(" ","")
                cadena = cadena.strip()

                if st.button('Búscame'):
                    st.write('---')
                    for indice,  (clave, valor) in enumerate(dic_query.items()):
                        st.toast(f'Consultando cromosoma {clave}', icon='🕠')
                        h = time.time()
                        seq = fst.fetch_genomic_dna_online(clave,valor[0], valor[1],st.session_state.email)
                        j = time.time()
                        st.success(f'Cromosoma {clave} obtenido en {round(j-h,2)} segundos. Desde posición {valor[0]} hasta {valor[1]}')
                        encontrados = fst.buscar_secuencia_en_fasta_todo(seq, cadena)
                        
                        if len(encontrados)>0:
                            encontrados = conversion_tuplas(dic_query, clave, encontrados)
                            with st.expander('Consulta obtenida del NCBI.'):
                                st.write('''Sólo se muestran 700 caracteres\n''', f'\n{seq[:700]}')
                            mostrar_resultados_acotados(indice, encontrados)
                        else:
                            st.error(f'Secuencia \'{cadena}\' no encontrada en cromosoma {clave}')
                        st.write('---')



        # Búsqueda OffLine
        if st.session_state.data_type =='OffLine':
            archivos_fasta = listar_archivos_fasta(ruta_crom_offline)
            nombres_sin_extension = [clave for clave, _ in archivos_fasta.items()]
            st.write('---')
            st.write('### Elige los cromosomas')
            seleccionados = st.multiselect('Elige Cromosomas', sorted(nombres_sin_extension), placeholder='Seleccionar Cromosa')
            
            if seleccionados:
                st.write('''##### Indica el inicio y final de la secuencia del cromosoma ''')
                dic_query = seleccion_rangos(seleccionados)
                st.write('''##### Indica la cadena que quieres buscar''')
    
                cadena = st.text_input('Introduce cadena a buscar').upper().replace(" ","")
                cadena = cadena.strip()
                
        
                if st.button('Busca'):
                    st.write('---')
                    for indice, (clave, valor) in enumerate(dic_query.items(), 0):
                        st.toast(f'Buscando cadena = {cadena} en cromosoma {clave}',icon='🕵️‍♂️')
                        encontrados = fst.buscar_secuencia_en_fasta_file_todo_acotado(archivos_fasta[clave][1],cadena, valor[0], valor[1])
                        
                        if len(encontrados)>0:
                            st.write(f'#### Resultados para {clave}')
                            mostrar_resultados_acotados(indice, encontrados)
                        else:
                            st.write(f'#### Resultados para {clave}')
                            st.error(f'Secuencia \'{cadena}\' no encontrada en cromosoma {clave}')
                        st.write('---')
            
#--------------------------------------SUBSTITUIR -------------------------------------------------------  
#--------------------------------------SUBSTITUIR -------------------------------------------------------
                     
    if page == 'Substituir':
            st.sidebar.write("""
            ## Sobre:
            Modificaciones de las secuencias con la enfermedad
                    """)
            
            st.subheader("Crea archivo con secuencia substituida")
            subcolom1, subcolom2 = st.columns(2)
            with subcolom1:
                 st.session_state.data_type = st.radio("Elige una opción:", ('Cadena',  'Cromosoma'), index = 1)
            
#--------------------------------------SUBSTITUIR -CADENA------------------------------------------------------            
            if st.session_state.data_type =='Cadena':
                st.write('---')
                st.write('#### Introduce Secuencia de Referencia')
                sequ = st.text_input('Introduce secuencia de referencia',placeholder='Pega o escribe aquí la secuencia').upper().replace(" ","")
                sequ = sequ.strip()
                
                if sequ:
                    st.write('#### Introduce cadena para buscar')
                    cadena = st.text_input('Introduce cadena',placeholder='Escribe la cadena que se buscará').upper().replace(" ","")
                    cadena = cadena.strip()
                    st.write('#### Indica el rango de la secuencia de referencia')
                    acotado = seleccion_rango_busqueda_cadena(len(sequ))
                    posiciones = []

                    st.write('---')
                    st.write('##### Ver resultados')

                    with st.expander('## Ver resultados obtenidos'):
                        if cadena != '':
                            posiciones = buscar_secuencia(sequ,cadena, acotado[0], acotado[1])
                            
                            if posiciones == []:
                                st.error('No encontrada')
                            else:
                                mostrar_resultados_acotados_EXPERIMENTAL(1,posiciones)
                            if 'result1' in st.session_state: 
                                if not st.session_state['result1'].empty:
                                    num_res = len(st.session_state['result1'])
                                    st.success(f'Resultados encontrados = {num_res}')
                                    st.dataframe(st.session_state['result1'], width=600)
                    
                    if cadena != '' and posiciones != []:
                        col1,col2 = st.columns(2)
                        
                        with col1:
                            st.write("#### Elige Opción")
                            st.session_state.data_accion = st.radio(" ", ('Eliminar','Substituir'), index=0)
                            if  st.session_state.data_accion=='Eliminar':
                                dic_ac = crea_diccionario_para_eliminar_cadena(posiciones)
                                st.write('---')
                                st.success(fst.secuencia_modificada(sequ, posiciones,dic_ac))
                        with col2:
                            if st.session_state.data_accion == 'Substituir':
                                st.write("#### Indica Cadena substituta")
                                subs = st.text_input('Introduce cadena substituta',placeholder='Escribe aquí la cadena substituta').upper().replace(" ","")
                                subs = subs.strip()
                                dic_ac = crea_diccionario_para_substituir_cadena(posiciones,subs)

                        if cadena != '' and posiciones != [] and st.session_state.data_accion == 'Substituir' and subs != '':
                            st.write('---')
                            st.write("#### Resultado de substituir Cadena")
                            st.success(f'{fst.secuencia_modificada(sequ, posiciones,dic_ac)}')



#--------------------------------------SUBSTITUIR -CROMOSOMA------------------------------------------------------
            if st.session_state.data_type =='Cromosoma':
                archivos_fasta = listar_archivos_fasta(ruta_crom_offline)
                #print(archivos_fasta)
                nombres_sin_extension = [clave for clave, _ in archivos_fasta.items()]
                st.write('---')
                st.write('### Elige Enfermedad')
                opcion_seleccionada = st.selectbox("Selecciona una enfermedad", sorted(list(indices_enfermedades_clndn.keys())), index=None)
                if opcion_seleccionada:
                    crom_afect = busca_crom_en_arch(opcion_seleccionada)
                    crom_afect = sorted(crom_afect, key=clave_ordenacion)

                    dr = ', '.join(c for c in crom_afect)
                    st.write(f'Cromosomas afectados por {opcion_seleccionada}: {dr}')

                    conj1 = fst.convertir_nombre_crom_del_vcf_a_ncbi(crom_afect)
                    conj2 = set(nombres_sin_extension)
                    interseccion = conj1.intersection(conj2)
                    lista1 = list(interseccion)
                    lista1.sort()

                    st.write('---')
                    st.write('### Elige cromosoma')
                    seleccionados = st.selectbox('Elige Cromosomas',lista1,index=None ,placeholder='Seleccionar Cromosa')
                    if seleccionados:
                        st.write('### Elige el rango')
                        acotado = seleccion_rango_cromosoma(seleccionados)

                        _, busqueda = busca_en_arch(opcion_seleccionada,diccionario_correspondencia.get(seleccionados),acotado[0],acotado[1])
                        st.toast(f'Búsqueda finalizada', icon='✅')

                        df = mostrar_result(busqueda,columnas_busqueda) 
                        if df.empty: 
                            st.error(f'La enfermedad {opcion_seleccionada} no se encuentra en el cromosoma {diccionario_correspondencia.get(seleccionados)}')
                        else:
                            st.write(f'### Posiciones del cromosoma {diccionario_correspondencia.get(seleccionados)} afectadas')
                            with st.expander(label='Haz clic para ver las posiciones afectadas'):
                                st.success(f'Número de resultados = {len(df)}')
                                st.dataframe(df,width=950)
                            # creo diccionario para tener todas las claves
                            dic_aux_pos_uni, dic_aux_pos_repe = crea_dic_pos_sub(busqueda)

                            claves_unidas = sorted(set(dic_aux_pos_uni.keys()) | set(dic_aux_pos_repe.keys()))

                            #print(claves_unidas)
                            lista_enteros = [int(elemento) for elemento in claves_unidas]
                            lista_enteros_ordenada = sorted(lista_enteros)
                            claves_unidas = [str(elemento) for elemento in lista_enteros_ordenada]
                                
                            #print(claves_unidas)
                            st.write("#### Elige tipo substitución")
                            st.write("Para seleccionar posiciones específicas, pulsa 'SELECCION'  y elige las posiciones deseadas de la lista a continuación.\nSi deseas seleccionar todas las posiciones, simplemente haz clic en el botón 'Todo'")
                                
                            col1,col2 = st.columns(2)
                            with col1:
                                st.session_state.data_sub = st.radio(" ", ('Seleccion','Todo'), index=0)

                            if st.session_state.data_sub == 'Todo':
                            # Multiselect para que el usuario elija las posiciones a cambiar
                                with st.expander(label='Muestra todas las posiciones seleccionadas'):
                                    selec = st.multiselect('Elige las posiciones que quieres cambiar', claves_unidas, placeholder='Selecciona Pos', default=claves_unidas)
                            
                            
                            if st.session_state.data_sub == 'Seleccion':
                                selec = st.multiselect('Elige las posiciones que quieres cambiar', claves_unidas, placeholder='Selecciona Pos',default=None)
                            st.write('---')
                                
                            proceso_seleccion(selec, dic_aux_pos_repe, dic_aux_pos_uni, ruta_crom_offline, seleccionados, acotado)
    
    #--------------------------------------kmer -CROMOSOMA------------------------------------------------------

    
    if page == 'K-mer':
        st.sidebar.write("""
            ## Sobre:
            k-mers \n
                    """)
        st.title('Extracción y análisis de kmers')
        archivos_fasta = listar_archivos_fasta(ruta_crom_offline)
               
        nombres_sin_extension = [clave for clave, _ in archivos_fasta.items()]
        subcolumna1, subcolumna2 = st.columns(2)
        
        with subcolumna1:
            st.session_state.data_type = st.radio("Elige una opción", ('k-mer','Grupo_Kmers', 'Comparativa'), index = 0)

        if st.session_state.data_type =='k-mer':
            seleccionados = st.selectbox('Elige Cromosomas', sorted(nombres_sin_extension),index=None ,placeholder='Seleccionar Cromosa')
            if seleccionados:
                acotado = seleccion_rango_cromosoma(seleccionados)
                valor_k = st.number_input("Valor de k", min_value = 1, value = 14)
                #st.success(f'{archivos_fasta.get(seleccionados)[1]}')
            if st.button('Ejecutar análisis'):
                st.write('---')
                sequ = fst.extraer_subsecuencia_de_fasta_file(archivos_fasta.get(seleccionados)[1],acotado[0],acotado[1])
                h = time.time()
        
                diccionario_kmers,diccionario_multiplicidades,kmers_agrupados_por_multiplicidad,frecuencias_relativas = fst.procesar_kmers(sequ, valor_k)
                j = time.time()
  

                st.success(f'Kmers Generados = {len(diccionario_kmers)}')
           

                # Filtrar el diccionario para incluir solo multiplicidades > 1
                diccionario_multiplicidades_filtrado = {multiplicidad: num_kmers 
                                            for multiplicidad, num_kmers in diccionario_multiplicidades.items() 
                                            if multiplicidad > 80}
                # Prepara datos para la gráfica
                multiplicidades = list(diccionario_multiplicidades_filtrado.keys())
                
                numero_de_kmers = list(diccionario_multiplicidades_filtrado.values())
                    # Crea figura y ejes para la gráfica
                fig, ax = plt.subplots(figsize=(40, 15))  
                ax.bar(multiplicidades, numero_de_kmers, color='skyblue')
                ax.set_xlabel('Multiplicidad (Número de veces que un k-mer aparece)',fontsize=54)
                ax.set_ylabel('Número de k-mers con esa multiplicidad',fontsize=54)
                ax.set_title('Distribución de multiplicidades de k-mers',fontsize=74)
                ax.set_xlim(-0.001, 400)  
                plt.xticks(rotation=45, fontsize=48)
                plt.yticks(fontsize=48)
                plt.tight_layout()

                # Mostrar la gráfica en Streamlit pasando la figura 
                st.write('##### Número de kmers por multiplicicdad')
                st.pyplot(fig)

                st.write('---')
                st.write('##### Frecuencias relativas')
                h = time.time()
                mostrar_grafico(frecuencias_relativas)
                j = time.time()
                print(f'Tiempo de ejecución = {round(j-h,2)}')

                entropias = calcular_entropia_cada_kmers(diccionario_kmers)

 
                st.write('---')
                st.write('##### Entropía para Cada k-mer ')
                texto_html = f"<p style='color: blue;'>Sirve para identificar k-mers que tienen características inusuales, como aquellos que son más o menos frecuentes de lo esperado, lo que podría indicar regiones funcionales específicas o mutaciones.</p>"
                st.markdown(texto_html, unsafe_allow_html=True)
                mostrar_grafico2(entropias,'Entropia')


                # Limitar a un máximo de 100 k-mers por cada clave
                kmers_agrupados_por_multiplicidad_limitado = {k: v[:50] for k, v in kmers_agrupados_por_multiplicidad.items()}

                # Crear los datos para el DataFrame
                data = {
                    "Número": list(kmers_agrupados_por_multiplicidad_limitado.keys()),
                    "Secuencias": list(kmers_agrupados_por_multiplicidad_limitado.values())
                }


                df = pd.DataFrame(data)

                st.write('---')
                st.write('##### Extracto de Secuencias de K-mers por multiplicidad')
                st.dataframe(df,width=700)
                st.write('---')

                hhh = calcular_entropia_kmers(diccionario_kmers)
                st.success(hhh)




        if st.session_state.data_type =='Grupo_Kmers':
            st.write('---')
            st.write('### Elige Enfermedad')
            opcion_seleccionada = st.selectbox("Selecciona una enfermedad", sorted(list(indices_enfermedades_clndn.keys())), index=None)
            if opcion_seleccionada:
                crom_afect = busca_crom_en_arch(opcion_seleccionada)
                crom_afect = sorted(crom_afect, key=clave_ordenacion)
                dr = ', '.join(c for c in crom_afect)
                st.write(f'Cromosomas afectados por {opcion_seleccionada}: {dr}')

                conj1 = fst.convertir_nombre_crom_del_vcf_a_ncbi(crom_afect)
                conj2 = set(nombres_sin_extension)
                interseccion = conj1.intersection(conj2)
                lista1 = list(interseccion)
                lista1.sort()

                st.write('---')
                st.write('### Elige cromosoma')
                seleccionados = st.selectbox('Elige Cromosomas',lista1,index=None ,placeholder='Seleccionar Cromosa')
                if seleccionados:
                    acotado = seleccion_rango_cromosoma(seleccionados)

                    _, busqueda = busca_en_arch(opcion_seleccionada,diccionario_correspondencia.get(seleccionados),acotado[0],acotado[1])
                    if seleccionados:
                        st.write(f'### Posiciones del cromosoma {diccionario_correspondencia.get(seleccionados)} afectadas')
                        df = mostrar_result(busqueda,columnas_busqueda) 
                        with st.expander(label='Haz clic para ver las posiciones afectadas'):
                            st.success(f'Número de resultados = {len(df)}')
                            st.dataframe(df,width=950)
                        st.write('---')
                        st.write('### Indica el número de nucleótidos de diferencia entre los grupos')

                        diferencia = st.number_input(label=f'Indica la diferencia', min_value=9000, max_value=1500000, value=None)
                        # creo diccionario para tener todas las claves
                        dic_aux_pos_uni, dic_aux_pos_repe = crea_dic_pos_sub(busqueda)

                        claves_unidas = sorted(set(dic_aux_pos_uni.keys()) | set(dic_aux_pos_repe.keys()))

                        #print(claves_unidas)
                        lista_enteros = [int(elemento) for elemento in claves_unidas]
                        lista_enteros_ordenada = sorted(lista_enteros)
                        claves_unidas = [str(elemento) for elemento in lista_enteros_ordenada]
                        if diferencia:
                            
                            lista_de_grupos = fst.distancia_entre_grupos(claves_unidas,diferencia)
                            


                            st.success(f'Número de grupos formados {len(lista_de_grupos)}')
                            mostrar_grupos(lista_de_grupos,seleccionados,dic_aux_pos_uni,dic_aux_pos_repe,archivos_fasta,seleccionados,acotado)
        
        if st.session_state.data_type =='Comparativa':
             # Ejemplo de uso del método
            directorio_actual = 'ANALISIS'  
            info_archivos = explorar_directorio_y_extraer_informacion(directorio_actual)
            archivos_fasta = listar_archivos_fasta(ruta_crom_offline)
            #print(archivos_fasta)
            if info_archivos:
                numeros_grupo = [dic['grupo'] for dic in info_archivos]
                numeros_grupo = sorted(numeros_grupo, key=int)
                st.write('---')
                st.write('#### Elige grupo para analizar')
                grupo = st.selectbox('Elige Grupo',numeros_grupo,index=None ,placeholder='Seleccionar grupo')
                if grupo:
                    for dic in info_archivos:
                        if dic['grupo'] == grupo:
                            st.write(f"###### Información del Grupo {grupo}:")
                            ini = dic['posicion_inicial']
                            fin = dic['posicion_final']
                            cromosoma = dic['nc_nombre']
                            #porcentaje = ((int(fin) - int(ini)+1) / int(acotado[1])) * 100
                            st.warning(f'Longitud de nucleótidos que se analizarán = {int(fin) - int(ini)+1} que corresponde al cromosoma {cromosoma}')
                            #print(f"  Ruta Completa: {dic['ruta_completa']}")
                            st.write(f"Posición Inicial: {ini}")
                            st.write(f"Posición Final: {fin}")
                            st.write(f"  Cromosoma: {dic['nc_nombre']}")


                    st.write(f" ### Tamaño del kmer")
                    k_val_max = st.number_input(label= f'Tamaño del K-mer',min_value=1, max_value=15,value=14,key='min')
                    k_val_min = st.number_input(label= 'Tama Min de K-mer',min_value=1, max_value=14,value=14,key='max')
                    #filtro = st.number_input(label= 'Filtro gráfica ',min_value=1, max_value=10,value=None)
                    # Configuración inicial
                    st.write('---')
                    st.write(f" ### Eliminar kmer no significativoS gráficas")
                    # Mostrar los valores seleccionados en inputs numéricos
                    min_value_input_fre,max_value_input_fre,min_value_input_ent,max_value_input_ent = intervalo_elimina_ruido()
                    st.write('---')




                    if st.button('Comparar')and k_val_min and k_val_max:
                        seq_mod = fst.leer_secuencia_de_fasta_file(dic['ruta_completa'])
                        seq_ori = fst.leer_secuencia_de_fasta_file_acotada(archivos_fasta.get(dic['nc_nombre'])[1],int(ini),int(fin))
                        for i in range(k_val_min,k_val_max+1):
                            st.write('---')
                            text= f'Generando estadísticas para K-mers de longitud = {i} en una secuencia de {int(fin) - int(ini)+1} nucleótidos'
                            st.markdown(f"""
                                        <div style='background-color: #e6f2ff; color: black; padding: 13px; border-radius: 6px;font-weight: bold;'>
                                            {text}
                                        </div>
                                        """, unsafe_allow_html=True)
                            #st.write('---')
                            #st.success(f'Mostrando Resultados de Cadena Mutada')
                            diccionario_kmers_mod,_,_,frecuencias_relativas_mod = fst.procesar_kmers(seq_mod, i)
                            #mostrar_grafico2(frecuencias_relativas_mod,'Frecuencia')
                            entropias_mod = calcular_entropia_cada_kmers(diccionario_kmers_mod)
                            #mostrar_grafico2(entropias_mod,'Entropia')

                            #st.write('---')
                            diccionario_kmers_ori,_,_,frecuencias_relativas_ori = fst.procesar_kmers(seq_ori, i)
                            #mostrar_grafico2(frecuencias_relativas_ori,'Frecuencia')
                            entropias_ori = calcular_entropia_cada_kmers(diccionario_kmers_ori)
                            #mostrar_grafico2(entropias_ori,'Entropia')
                            st.write('---')
                            st.success(f'Comparativa de frecuencia relativa de k-mers')
                            comp_freq_kmers = fst.comparar_diccionarios(frecuencias_relativas_ori,frecuencias_relativas_mod)
                            comp_entro_kmers = fst.comparar_diccionarios(entropias_ori,entropias_mod)
                            df_fre = mostrar_grafico3(comp_freq_kmers,'Frecuencia_Comp',min_value_input_fre,max_value_input_fre)
                            st.write('---')
                            st.success(f'Comparativa de entropias relativa de k-mers')
                            def_entr = mostrar_grafico3(comp_entro_kmers,'Entropia_Comp',min_value_input_ent,max_value_input_ent)
                            def_entr.to_csv(f'Output/kmers_entorpias g{grupo}_k{i}_{cromosoma}_{ini}_{fin}.csv', index=False)
                            df_fre.to_csv(f'Output/kmers_frec_g{grupo}_k{i}_{cromosoma}_{ini}_{fin}.csv', index=False)
                            st.success(f'Archivo den entropias guardado en Output/kmers_entorpias_g{grupo}_k{i}_{cromosoma}_{ini}_{fin}.csv')
                            st.success(f'Archivo den frec guardado en Output/kmers_frec_g{grupo}_k{i}_{cromosoma}_{ini}_{fin}.csv')
                            
                            df_fre.sort_values(by='Frecuencia_Comp', ascending=False, inplace=True)
                            
                            top_100_kmers = df_fre.head(50)

                            # Título de la aplicación
                            st.markdown('###### Top 50 k-mers por Frecuencia Comparativa Mayor a menor')
                            with st.expander(label='Top 50 k-mers por Frecuencia Comparativa Mayor a menor'):
                                # Mostrar el DataFrame en Streamlit
                                st.dataframe(top_100_kmers,width=800)

                            df_fre.sort_values(by='Frecuencia_Comp', ascending=True, inplace=True)
                            top_100_kmers = df_fre.head(50)

                            st.markdown('###### Top 50 k-mers por Frecuencia Comparativa menor a Mayor')
                            with st.expander(label='Top 50 k-mers por Frecuencia Comparativa menor a Mayor'):
                                # Mostrar el DataFrame en Streamlit
                                st.dataframe(top_100_kmers,width=800)

                            #df_fre.sort_values(by='Entropia_Comp', ascending=False, inplace=True)
                            ###########################################
                            def_entr.sort_values(by='Entropia_Comp', ascending=False, inplace=True)
                            
                            top_100_kmers = def_entr.head(50)

                            # Título de la aplicación
                            st.markdown('###### Top 50 k-mers por entropía Comparativa Mayor a menor')

                            # Mostrar el DataFrame en Streamlit
                            with st.expander(label='Top 50 k-mers por entropía Comparativa mayor a menor'):
                                st.dataframe(top_100_kmers,width=800)

                            def_entr.sort_values(by='Entropia_Comp', ascending=True, inplace=True)
                            top_100_kmers = def_entr.head(50)
                            st.markdown('###### Top 50 k-mers por entropía Comparativa menor a Mayor')
                            with st.expander(label='Top 50 k-mers por entropía Comparativa menor a mayor'):
                            # Mostrar el DataFrame en Streamlit
                                st.dataframe(top_100_kmers,width=800)

                        #st.warning('asdfasdfasdf')
@st.cache_data                        
def mostrar_analis_kmer_top(df_fre,def_entr):
    df_fre.sort_values(by='Frecuencia_Comp', ascending=False, inplace=True)
                            
    top_100_kmers = df_fre.head(50)

                            # Título de la aplicación
    st.markdown('###### Top 50 k-mers por Frecuencia Comparativa Mayor a menor')
    with st.expander(label='Top 50 k-mers por Frecuencia Comparativa Mayor a menor'):
        # Mostrar el DataFrame en Streamlit
        st.dataframe(top_100_kmers,width=800)

    df_fre.sort_values(by='Frecuencia_Comp', ascending=True, inplace=True)
    top_100_kmers = df_fre.head(50)
    st.markdown('###### Top 50 k-mers por Frecuencia Comparativa menor a Mayor')
    with st.expander(label='Top 50 k-mers por Frecuencia Comparativa menor a Mayor'):
                            # Mostrar el DataFrame en Streamlit
        st.dataframe(top_100_kmers,width=800)
                        #df_fre.sort_values(by='Entropia_Comp', ascending=False, inplace=True)
                        ###########################################
    def_entr.sort_values(by='Entropia_Comp', ascending=False, inplace=True)
                            
    top_100_kmers = def_entr.head(50)

                            # Título de la aplicación
    st.markdown('###### Top 50 k-mers por entropía Comparativa Mayor a menor')

                            # Mostrar el DataFrame en Streamlit
    with st.expander(label='Top 50 k-mers por entropía Comparativa mayor a menor'):
        st.dataframe(top_100_kmers,width=800)

    def_entr.sort_values(by='Entropia_Comp', ascending=True, inplace=True)
    top_100_kmers = def_entr.head(50)
    st.markdown('###### Top 50 k-mers por entropía Comparativa menor a Mayor')
    with st.expander(label='Top 50 k-mers por entropía Comparativa menor a mayor'):
                            # Mostrar el DataFrame en Streamlit
        st.dataframe(top_100_kmers,width=800)
def mostrar_grupos(lista_de_grupos,cromosoma,dic_aux_pos_uni,dic_aux_pos_repe,archivos_fasta,seleccionados,acotado):
    dic_grupo_informacion ={}

    for indice, grupo in enumerate(lista_de_grupos,1):
        text= f'Elige posiciones para el grupo {indice}'
        st.markdown(f"""
            <div style='background-color: #e6f2ff; color: black; padding: 10px; border-radius: 5px;font-weight: bold;'>
                {text}
            </div>
            """, unsafe_allow_html=True)
        #st.write(f'##### Elige posiciones para el grupo {indice}')
        with st.expander(label='Muestra todas las posiciones seleccionadas'):
            selec = st.multiselect('Elige las posiciones que quieres cambiar', grupo, placeholder='Selecciona Pos', default=grupo,key=indice)
           
            
            rango = seleccion_rango_grupos(indice,selec,cromosoma)
            
            
            if selec:
                    # Filtrar los diccionarios basados en la selección y comprobar las posiciones repetidas y únicas
                    dic_selec_repe = {item: dic_aux_pos_repe[item] for item in selec if item in dic_aux_pos_repe}
                    dic_selec_uni = {item: dic_aux_pos_uni[item] for item in selec if item in dic_aux_pos_uni}

                    # Inicializar diccionario para posiciones con múltiples alternativas
                    dic_selec_repe_uni = {}
                    if dic_selec_repe:
                        st.error('#### Existen posiciones que tienen varias alternativas')
                        # Procesar posiciones con más de una alternativa
                        dic_selec_repe_uni = posiciones_mas_alternativas(dic_selec_repe)

                    
                    # Unir los diccionarios de selecciones únicas y repetidas
                    diccionario_unido = {**dic_selec_uni, **dic_selec_repe_uni}


                    # Ordenar las claves del diccionario unido y crear un nuevo diccionario ordenado
                    #nuevo = {clave: diccionario_unido[clave] for clave in sorted(diccionario_unido)}

                    # Ordenar las claves convirtiéndolas a enteros antes de compararlas
                    claves_ordenadas = sorted(diccionario_unido, key=lambda x: int(x))

                    # Crear un nuevo diccionario ordenado
                    diccionario_unido = {clave: diccionario_unido[clave] for clave in claves_ordenadas}

            

                    # Buscar solapamientos entre las posiciones seleccionadas
                    dic_solapes = fst.busca_solapamientos_posiciones(diccionario_unido)


                    #print(dic_solapes)
                    
                    #ruta_archivo_modificado = f'Output/{seleccionados}_mod.fasta'
                    
                    # Verificar y mostrar resultados de solapamientos
                    dic_solapes, diccionario_unido = verifica_y_muestra_tab(dic_solapes, diccionario_unido)
                    dic_grupo_informacion[str(indice)]= {'dic_unido':diccionario_unido,'posiciones': rango }

    if dic_grupo_informacion:
        clave_grupo = list(dic_grupo_informacion.keys())
        #print(clave_grupo)
        st.write('---')
        grupos_seleccionados = st.multiselect('Elige los grupos a crear', clave_grupo, placeholder='Selecciona Grupo', default=clave_grupo)

        if st.button('Modficiar y Guardar'):
            rutas_para_archivos = crear_estructura_directorios_y_subcarpetas(grupos_seleccionados)
            #print(rutas_para_archivos)
            for claves, diccs in dic_grupo_informacion.items():
                # Preparar el diccionario de sustituciones
                dic_subs = diccs['dic_unido']
                dic_subs = {
                (int(clave), int(clave) + len(valor[0][0]) - 1): (valor[0][0], valor[0][1])
                for clave, valor in dic_subs.items()
                }
                print('afadfadsadfadsfasdfadadfadadfasdfdafdadfsasdfasdf')
                #print(dic_subs)


                sequen = fst.modificar_fasta_solo(archivos_fasta.get(seleccionados)[1],dic_subs)
                
                nom_arc = f'{rutas_para_archivos[claves]}/{seleccionados}_{claves}'
                acotamiento = diccs['posiciones']
                sequen= fst.extraer_subsecuencia_de_fasta_online(sequen,acotamiento[0],acotamiento[1])
                guardar_archivo_fasta_kmer(nom_arc, acotamiento[0], acotamiento[1],sequen)
            st.session_state.guardado = True
        if 'guardado' in st.session_state:
             
            # Ejemplo de uso del método
            directorio_actual = 'ANALISIS'  
            info_archivos = explorar_directorio_y_extraer_informacion(directorio_actual)

            # Imprimir la información extraída para verificar

            
            if info_archivos:
                numeros_grupo = [dic['grupo'] for dic in info_archivos]
                numeros_grupo = sorted(numeros_grupo, key=int)
                st.write('---')
                st.write('#### Elige grupo para analizar')
                grupo = st.selectbox('Elige Grupo',numeros_grupo,index=None ,placeholder='Seleccionar grupo')
                
                if grupo:
                    for dic in info_archivos:
                        if dic['grupo'] == grupo:
                            st.write(f"###### Información del Grupo {grupo}:")
                            ini = dic['posicion_inicial']
                            fin = dic['posicion_final']
                            porcentaje = ((int(fin) - int(ini)+1) / int(acotado[1])) * 100
                            st.warning(f'Longitud de nucleótidos que se analizarán = {int(fin) - int(ini)+1} que corresponde a {porcentaje:.2f} % del cromosoma {cromosoma}')
                            #print(f"  Ruta Completa: {dic['ruta_completa']}")
                            st.write(f"Posición Inicial: {ini}")
                            st.write(f"Posición Final: {fin}")
                            st.write(f"  Cromosoma: {dic['nc_nombre']}")
                            
                    st.write(f" ### Tamaño máximo del kmer")
                    k_val_max = st.number_input(label= f'Tama Max de K-mer',min_value=1, max_value=15,value=14,key='min')
                    st.write(f" ### Tamaño min del kmer")
                            
                    k_val_min = st.number_input(label= 'Tama Min de K-mer',min_value=1, max_value=14,value=14,key='max')
                    st.write(f" ### Tamaño ventana Exploración")
                    tam_ventana = st.number_input(label= 'Tama ventana exploración',min_value=1, max_value=100000,value=None)


                    filtro = st.number_input(label= 'Filtro gráfica ',min_value=1, max_value=10,value=None)
                    
                    if filtro:
                        seq = fst.leer_secuencia_de_fasta_file(dic['ruta_completa'])
                        if tam_ventana:
                            if st.button('Analizando'):
                                longitud_secuencia = len(seq)
                                paso = tam_ventana
                                aux_pos = 0
                                for i in range(0, longitud_secuencia, paso): 
                                    longitud_sub_seq = len(seq) - aux_pos + 1
                                    for m in range(k_val_min,k_val_max+1):
                                        st.write('---')
                                        text= f'Generada estadísticas para K-mers de longitud = {m} en una secuencia de {longitud_sub_seq} nucleótidos'
                                        st.markdown(f"""
                                    <div style='background-color: #e6f2ff; color: black; padding: 13px; border-radius: 6px;font-weight: bold;'>
                                        {text}
                                    </div>
                                    """, unsafe_allow_html=True)
                                        st.write('---')
                                        h = time.time()
                                        diccionario_kmers,diccionario_multiplicidades,kmers_agrupados_por_multiplicidad,frecuencias_relativas = fst.procesar_kmers(seq, i)
                                        j = time.time()
                                        
                                        #diccionario_multiplicidades = fst.contar_multiplicidades(diccionario_kmers)
                                        # Filtrar el diccionario para incluir solo multiplicidades > 1
                                        diccionario_multiplicidades_filtrado = {multiplicidad: num_kmers 
                                                                    for multiplicidad, num_kmers in diccionario_multiplicidades.items() 
                                                                        if multiplicidad > filtro}
                                            # Prepara datos para la gráfica
                                        multiplicidades = list(diccionario_multiplicidades_filtrado.keys())
                                            
                                        numero_de_kmers = list(diccionario_multiplicidades_filtrado.values())
                                                # Crea figura y ejes para la gráfica
                                        fig, ax = plt.subplots(figsize=(40, 15))  
                                        ax.bar(multiplicidades, numero_de_kmers, color='skyblue')
                                        ax.set_xlabel('Multiplicidad (Número de veces que un k-mer aparece)',fontsize=54)
                                        ax.set_ylabel('Número de k-mers con esa multiplicidad',fontsize=54)
                                        ax.set_title('Distribución de multiplicidades de k-mers',fontsize=74)
                                        ax.set_xlim(-0.001, 400)  
                                        plt.xticks(rotation=45, fontsize=48)
                                        plt.yticks(fontsize=48)
                                        plt.tight_layout()
                                        st.pyplot(fig)
                                        plt.close(fig)

                                        h = time.time()
                                        mostrar_grafico2(frecuencias_relativas,'Frecuencia')
                                        j = time.time()
                                        print(f'Tiempo de ejecución = {round(j-h,2)}')
                                        entropias = calcular_entropia_cada_kmers(diccionario_kmers)

                                        h = time.time()
                                        mostrar_grafico2(entropias,'Entropia')
                                        j = time.time()
                                        print(f'Tiempo de ejecución = {round(j-h,2)}')
                                    aux_pos+= paso
                        else:
                            if st.button('Analizando'):
                                for i in range(k_val_min,k_val_max+1):
                                    h = time.time()
                                    diccionario_kmers,diccionario_multiplicidades,kmers_agrupados_por_multiplicidad,frecuencias_relativas = fst.procesar_kmers(seq, i)
                                    j = time.time()
                                    print(f'Tiempo de ejecución =  {round(j-h,2)} segundos')
                                    st.write('---')
                                    text= f'Generada estadísticas para K-mers de longitud = {i}'
                                    st.markdown(f"""
                                    <div style='background-color: #e6f2ff; color: black; padding: 13px; border-radius: 6px;font-weight: bold;'>
                                        {text}
                                    </div>
                                    """, unsafe_allow_html=True)
                                    st.write('---')
                                    #st.success(f'Generado diccionar kmers con k = {i}')
                                    #diccionario_multiplicidades = fst.contar_multiplicidades(diccionario_kmers)

                                    # Filtrar el diccionario para incluir solo multiplicidades > 1
                                    diccionario_multiplicidades_filtrado = {multiplicidad: num_kmers 
                                                                for multiplicidad, num_kmers in diccionario_multiplicidades.items() 
                                                                if multiplicidad > filtro}
                                    # Prepara datos para la gráfica
                                    multiplicidades = list(diccionario_multiplicidades_filtrado.keys())
                                    
                                    numero_de_kmers = list(diccionario_multiplicidades_filtrado.values())
                                        # Crea figura y ejes para la gráfica
                                    fig, ax = plt.subplots(figsize=(40, 15))  
                                    ax.bar(multiplicidades, numero_de_kmers, color='skyblue')
                                    ax.set_xlabel('Multiplicidad (Número de veces que un k-mer aparece)',fontsize=54)
                                    ax.set_ylabel('Número de k-mers con esa multiplicidad',fontsize=54)
                                    ax.set_title('Distribución de multiplicidades de k-mers',fontsize=74)
                                    ax.set_xlim(-0.001, 400)  
                                    plt.xticks(rotation=45, fontsize=48)
                                    plt.yticks(fontsize=48)
                                    plt.tight_layout()
                                    st.pyplot(fig)

                                    h = time.time()
                                    mostrar_grafico2(frecuencias_relativas,'Frecuencia')
                                    j = time.time()
                                    print(f'Tiempo de ejecución = {round(j-h,2)}')

                                    entropias = calcular_entropia_cada_kmers(diccionario_kmers)

                                    h = time.time()
                                    mostrar_grafico2(entropias,'Entropia')
                                    j = time.time()
                                    print(f'Tiempo de ejecución = {round(j-h,2)}')


    
            
def calcular_entropia_cada_kmers(diccionario_kmers):
    # Calcula el total de k-mers observados
    total_kmers = sum(diccionario_kmers.values())

    # Diccionario para almacenar la entropía de cada k-mer
    entropia_kmers = {}

    # Iterar sobre cada k-mer y calcular su entropía
    for kmer, frecuencia in diccionario_kmers.items():
        # Calcular la frecuencia relativa de este k-mer
        p_x = frecuencia / total_kmers
        # Calcular la contribución de este k-mer a la entropía
        entropia_kmers[kmer] = -p_x * math.log2(p_x)

    return entropia_kmers

def calcular_entropia_kmers(diccionario_kmers):

    # Calcula la entropía de Shannon para el conjunto completo de k-mers:
    total_kmers = len(diccionario_kmers)
    frecuencias_relativas = {kmer: frec / total_kmers for kmer, frec in diccionario_kmers.items()}

    # Paso 3: Aplicar la fórmula de entropía de Shannon
    entropia = -sum(frecuencia * math.log2(frecuencia) for frecuencia in frecuencias_relativas.values())

    return entropia



def mostrar_grafico(frecuencias_relativas):


    df = pd.DataFrame(list(frecuencias_relativas.items()), columns=['kmer', 'Frecuencia'])

    # Filtrar el DataFrame para incluir solo los k-mers con frecuencia mayor que 1
    df = df[df['Frecuencia']>0.00005]

    # Ordenar el DataFrame por frecuencia para una visualización más clara

    df = df.sort_values(by='Frecuencia', ascending=False)

    # Crear un gráfico de barras con Altair
    chart = alt.Chart(df).mark_bar(size=2).encode(
        x='kmer:N',  # N indica categoría (nominal)
        y='Frecuencia:Q',  # Q indica cantidad (cuantitativa)
        tooltip=['kmer:N', 'Frecuencia:Q']  # Mostrar información al pasar el ratón
    ).properties(
        width=750,  # Ancho del gráfico
        height=400,  # Altura del gráfico
        title='Frecuencias de k-mers'  # Título del gráfico
    ).interactive()

    # Mostrar el gráfico en Streamlit
    st.altair_chart(chart, use_container_width=True)

    
def mostrar_grafico2(frecuencias_relativas, frecuenc):
    # Crear DataFrame con la columna de frecuencia ajustable
    df = pd.DataFrame(list(frecuencias_relativas.items()), columns=['kmer', frecuenc])

    # Filtrar el DataFrame para incluir solo los k-mers con una frecuencia por encima de un umbral
    df = df[df[frecuenc] > 0.00001]

    # Ordenar el DataFrame por la columna de frecuencia dinámica
    df = df.sort_values(by=frecuenc, ascending=False)


    # Crear un gráfico de barras con Altair usando la columna de frecuencia dinámica
    chart = alt.Chart(df).mark_bar(size=2).encode(
        x='kmer:N',  # N indica categoría (nominal)
        y=alt.Y(frecuenc + ':Q'),  # Q indica cantidad (cuantitativa), con frecuenc como columna dinámica
        tooltip=['kmer:N', alt.Tooltip(frecuenc + ':Q')]  # Mostrar información al pasar el ratón
    ).properties(
        width=750,  # Ancho del gráfico
        height=400,  # Altura del gráfico
        title=f'Frecuencias de k-mers usando {frecuenc}'  # Título del gráfico
    ).interactive()

    # Mostrar el gráfico en Streamlit
    st.altair_chart(chart, use_container_width=True)





def mostrar_grafico3(frecuencias_relativas, frecuenc,min_value_input_fre,max_value_input_fre):

    # Crear DataFrame con la columna de frecuencia ajustable
    df = pd.DataFrame(list(frecuencias_relativas.items()), columns=['kmer', frecuenc])

    # Tiempo de inicio para el filtrado

    min = "{:.6f}".format(min_value_input_fre)
    max = "{:.6f}".format(max_value_input_fre)

    # Filtrar el DataFrame para incluir solo los k-mers con una frecuencia por encima de un umbral
    df = df[~((df[frecuenc] >= float(min)) & (df[frecuenc] <= float(max)))]


    # Crear un gráfico de barras con Altair usando la columna de frecuencia dinámica
    chart = alt.Chart(df).mark_bar(size=2).encode(
        x='kmer:N',  # N indica categoría (nominal)
        y=alt.Y(frecuenc + ':Q'),  # Q indica cantidad (cuantitativa), con frecuenc como columna dinámica
        tooltip=['kmer:N', alt.Tooltip(frecuenc + ':Q')]  # Mostrar información al pasar el ratón
    ).properties(
        width=750,  # Ancho del gráfico
        height=400,  # Altura del gráfico
        title=f'Frecuencias de k-mers usando {frecuenc}'  # Título del gráfico
    ).interactive()

    # Mostrar el gráfico en Streamlit
    st.altair_chart(chart, use_container_width=True)
    return df

    



@st.cache_data
def agrupar_kmers_por_multiplicidad_cached(diccionario_kmers):
    return fst.agrupar_kmers_por_multiplicidad(diccionario_kmers)                 

def proceso_seleccion(selec, dic_aux_pos_repe, dic_aux_pos_uni, ruta_crom_offline, seleccionados, acotado):
    # Verificar si hay elementos seleccionados
    guardado = False
    if selec:
        # Filtrar los diccionarios basados en la selección y comprobar las posiciones repetidas y únicas
        dic_selec_repe = {item: dic_aux_pos_repe[item] for item in selec if item in dic_aux_pos_repe}
        dic_selec_uni = {item: dic_aux_pos_uni[item] for item in selec if item in dic_aux_pos_uni}

        # Inicializar diccionario para posiciones con múltiples alternativas
        dic_selec_repe_uni = {}
        if dic_selec_repe:
            st.error('#### Existen posiciones que tienen varias alternativas')
            # Procesar posiciones con más de una alternativa
            dic_selec_repe_uni = posiciones_mas_alternativas(dic_selec_repe)

        
        # Unir los diccionarios de selecciones únicas y repetidas
        diccionario_unido = {**dic_selec_uni, **dic_selec_repe_uni}


        # Ordenar las claves convirtiéndolas a enteros antes de compararlas
        claves_ordenadas = sorted(diccionario_unido, key=lambda x: int(x))

        # Crear un nuevo diccionario ordenado
        diccionario_unido = {clave: diccionario_unido[clave] for clave in claves_ordenadas}

        # Buscar solapamientos entre las posiciones seleccionadas
        dic_solapes = fst.busca_solapamientos_posiciones(diccionario_unido)

        st.write('#### Verificando si las posiciones seleccionadas se solapan')
        ruta_archivo_modificado = f'Output/{seleccionados}_mod.fasta'
        
        # Verificar y mostrar resultados de solapamientos
        dic_solapes, diccionario_unido = verifica_y_muestra(dic_solapes, diccionario_unido)

        
        if not dic_solapes:
            st.success(f'Ahora ya no hay solapamiento. Número de substituciones en total {len(diccionario_unido)}')
           

        # Opción para realizar las sustituciones
        if st.button('Substituir'):
            # Preparar el diccionario de sustituciones
            dic_subs = {
                (int(clave), int(clave) + len(valor[0][0]) - 1): (valor[0][0], valor[0][1])
                for clave, valor in diccionario_unido.items()
            }
            # Ordenar el diccionario por el primer elemento de la clave (tupla)
            claves_ordenadas = sorted(dic_subs.keys(), key=lambda x: x[0])

            # Crear un nuevo diccionario ordenado
            dic_subs = {clave: dic_subs[clave] for clave in claves_ordenadas}

            if st.session_state.data_sub == 'Seleccion':
                st.success(dic_subs)

            # Listar archivos FASTA disponibles
            archivos_fasta = listar_archivos_fasta(ruta_crom_offline)
            ruta_archivo_modificado = f'Output/{seleccionados}_mod.fasta'

            fst.modificar_fasta(archivos_fasta.get(seleccionados)[1], ruta_archivo_modificado, dic_subs)
            guardado = True


            with st.expander(label='Haz clic para ver las posiciones modificadas'):
                # Mostrar los cambios realizados
                for a, valor in dic_subs.items():
                    st.success(f'Cambio en la posición {a} cuyo valor anterior era {valor[0]}. Ahora es {valor[1]}')


        # Opción para descargar la secuencia FASTA modificada
        if 'mod' in st.session_state and st.session_state.data_sub !='Todo':
            st.download_button(
                label="Descargar Secuencia FASTA modificada",
                data=st.session_state['mod'],
                file_name="secuencia.fasta",
                mime="text/x-fasta"
            )
            # Mostrar las secuencias original y modificada
            if st.session_state.data_sub != 'Todo':
                ii = st.session_state['ori']
                ij = st.session_state['mod']
                st.write(f'{ii}')
                st.write(f'{ij}')
        else:
            st.write(f'Se guardará en {ruta_archivo_modificado}')
            if guardado:
                st.success('Se ha guardado ya')


def verifica_y_muestra_tab(dic_solapes, diccionario_unido):
    if dic_solapes:
        
        st.warning('##### ¡Existen solapamientos!')
        solap_selec, sola_no_selec = [],[]

        # Crear pestañas para cada solapamiento
        tabs = st.tabs([f'Solapamiento {i}' for i in range(1, len(dic_solapes) + 1)])

        for i, (clave, valor) in enumerate(dic_solapes.items()):
            with tabs[i]:  # Usar el tab correspondiente
                seleccion = st.selectbox('Elige una posición', valor, index=0, placeholder='Selecciona Pos', key=clave)
                solap_selec.append(seleccion)

                # Encontrar los no seleccionados eficientemente
                no_seleccionados = set(valor) - {seleccion}
                sola_no_selec.extend(no_seleccionados)

                # Mostrar información sobre los solapamientos
                for item in valor:
                    ref = diccionario_unido.get(str(item))[0][0]
                    alt = diccionario_unido.get(str(item))[0][1]
                    texto_html = f"<p style='color: blue;'>La posición {item} tiene como REF = {ref} que será subtituida por ALT = {alt}</p>"
                    st.markdown(texto_html, unsafe_allow_html=True)

        # Eliminar las claves no seleccionadas del diccionario principal
        for c in sola_no_selec:
            diccionario_unido.pop(str(c), None)  # Usar pop con None para evitar KeyError

        # Ordenar y actualizar el diccionario
        claves_ordenadas = sorted(diccionario_unido, key=lambda x: int(x))

        nuevo = {cada: diccionario_unido[cada] for cada in claves_ordenadas}
        dic_solapes = fst.busca_solapamientos_posiciones(nuevo)  # Asegúrate de que `fst` está definido/importado correctamente
      
    else:
        st.success('##### ¡No existen solapamientos!')
    return dic_solapes, diccionario_unido

def verifica_y_muestra(dic_solapes, diccionario_unido):
    if dic_solapes:
        
        st.warning('##### ¡Existen solapamientos!')
        solap_selec, sola_no_selec = [],[]

        
        # Procesar cada solapamiento y mostrar en expanders
        for indice, (clave, valor) in enumerate (dic_solapes.items(),1):
            with st.expander(label=f'Solapamiento {indice}'):
                # Registrar la selección del usuario
                seleccion = st.selectbox('Elige una posición', valor, index=0, placeholder='Selecciona Pos', key=clave)
                solap_selec.append(seleccion)
                # Encontrar los no seleccionados eficientemente
                no_seleccionados = set(valor) - {seleccion}
                sola_no_selec.extend(no_seleccionados)
                # Mostrar información sobre los solapamientos
                for item in valor:
                    ref = diccionario_unido.get(str(item))[0][0]
                    alt = diccionario_unido.get(str(item))[0][1]
                    texto_html = f"<p style='color: blue;'>La posición {item} tiene como REF = {ref} que será subtituida por ALT = {alt}</p>"
                    st.markdown(texto_html, unsafe_allow_html=True)
        # Eliminar las claves no seleccionadas del diccionario principal
        for c in sola_no_selec:
            diccionario_unido.pop(str(c), None)  # Usar pop con None para evitar KeyError

        # Ordenar y actualizar el diccionario
        claves_ordenadas = sorted(diccionario_unido,key=lambda x: int(x))

        nuevo = {cada: diccionario_unido[cada] for cada in claves_ordenadas}
        dic_solapes = fst.busca_solapamientos_posiciones(nuevo)
        
    return dic_solapes,diccionario_unido


def posiciones_mas_alternativas(dic_selec_repe):
    dic_repe_uni={}
    for clave, valor in dic_selec_repe.items():
        st.write(f'##### Posición {clave} tiene más de una alternativa. Elige una.')
        # Crear un selectbox con un estado que persista entre refrescos de la página
        dic_repe_uni[clave] = [st.selectbox(f'Elige una opción para {clave}. Significado (REF, ALT)',
                                                        dic_selec_repe[clave],
                                                        key=clave)]
        st.write('---')
    return dic_repe_uni

def mostrar_result(busqueda, columnas_m):
    st.toast(f'Mostrando Datos',icon='🧑‍💻')
    data = StringIO(busqueda)
    df = pd.read_csv(data, sep="\t", names=columnas_m, header=None, dtype={0: str})  
    #st.success(f'Número de resultados = {len(df)}')
    #st.dataframe(df,width=650)
    return df

@st.cache_data               
def crea_dic_pos_sub(busqueda):
    patron_pos1 = re.compile(r'^(?!#)[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)')
    dic_aux_pos_uni = {}
    dic_aux_pos_repe = {}
    for linea in busqueda.split('\n'):
        match = patron_pos1.match(linea)
        if match:
            if match.group(1) not in dic_aux_pos_uni:
                # caso de que exista diferentes alternativas separadas por ,
                if len(match.group(3).split(',')) > 1:
                    for cada in match.group(3).split(','):
                        dic_aux_pos_repe.setdefault(match.group(1), []).append((match.group(2), cada))
                else:
                    dic_aux_pos_uni[match.group(1)] = [(match.group(2), match.group(3) )]
            else:

                valor = dic_aux_pos_uni.pop([match.group(1)][0], None)[0]

                if valor is not None:
                    dic_aux_pos_repe.setdefault(match.group(1), []).append(valor)
                    dic_aux_pos_repe.setdefault(match.group(1), []).append((match.group(2), match.group(3)))
                else:

                    dic_aux_pos_repe.setdefault(match.group(1), []).append((match.group(2), match.group(3)))
    return dic_aux_pos_uni, dic_aux_pos_repe

@st.cache_data
def buscar_secuencia(sequ,cadena,ini,fin):
    return fst.buscar_secuencia_en_fasta_acotado(sequ,cadena,ini,fin)
       

def conversion_tuplas(dic_query, clave, encontrados):
    suma = dic_query[clave][0]-1
    #st.toast(dic_query[clave])
    aux_encon = []
    for elemento in encontrados:

        aux_encon.append((elemento[0]+suma, elemento[1]+suma))

    return aux_encon

def crea_diccionario_para_substituir_cadena(posiciones:list, cadena:str)->dict:
    dic_acciones = {}
    for cada in posiciones:
        dic_acciones[cada] = (f'{cadena}')
    return dic_acciones


@st.cache_data
def crea_diccionario_para_eliminar_cadena(posiciones:list)->dict:
    dic_acciones = {}
    for cada in posiciones:
        dic_acciones[cada] = ('')
    return dic_acciones

def mostrar_resultados_acotados_EXPERIMENTAL(indice, encontrados):
    aux = f'result{indice}'

    #st.success(f'{len(encontrados)} posiciones encontradas con la cadena {cadena} en el cromosoma {clave} acotado entre {valor[0]} y {valor[1]}')
    st.toast(f'Búsqueda finalizada. Creando DataFrame',icon='✅')
    columnas_encon=[ 'Inicio', 'Final']
    df = pd.DataFrame(encontrados,columns=columnas_encon)
    st.write('###### Resultados')
    st.session_state[aux] = df
    #st.dataframe(df, width=600)  
    


@st.cache_data
def mostrar_resultados_acotados_online(indice, encontrados, cadena, clave, valor):
    aux = f'result{indice}'

    st.success(f'{len(encontrados)} posiciones encontradas con la cadena {cadena} en el cromosoma {clave} acotado entre {valor[0]} y {valor[1]}')
    st.toast(f'Búsqueda finalizada. Creando DataFrame',icon='✅')
    columnas_encon=[ 'Cadena','Inicio', 'Final', 'Secuencia']
    df = pd.DataFrame(encontrados,columns=columnas_encon)

    st.dataframe(df, width=600)  
    st.session_state[aux] = df
    

def mostrar_resultados_acotados(indice, encontrados):
    aux = f'result{indice}'

    #st.success(f'{len(encontrados)} posiciones encontradas con la cadena {cadena} en el cromosoma {clave} acotado entre {valor[0]} y {valor[1]}')
    st.toast(f'Búsqueda finalizada. Creando DataFrame',icon='✅')
    columnas_encon=[ 'Inicio', 'Final']
    df = pd.DataFrame(encontrados,columns=columnas_encon)
    st.write('###### Resultados')
    st.dataframe(df, width=600)  
    st.session_state[aux] = df




def anyadir_fila_seq(df, seq, inicio, fin):
    nueva_fila = pd.DataFrame([[seq, inicio, fin]], columns=df.columns)
    return pd.concat([df, nueva_fila], ignore_index=True)

def anyadir_fila(df, cromosoma, inicio, fin, secuencia):
    # Añadir una fila al DataFrame
    nueva_fila = pd.DataFrame([[cromosoma, inicio, fin, secuencia]], columns=df.columns)
    #df = df.concat([df, nueva_fila], ignore_index=True)
    
    return pd.concat([df, nueva_fila], ignore_index=True)

         

def guardar_archivo_fasta(nombre, inicio, fin, seq):
    ruta = f'Output/{nombre}_{inicio}_{fin}.fasta'
    with open(ruta, 'w') as file:
        file.write(seq)

def guardar_archivo_fasta_kmer(nombre, inicio, fin, seq):
    ruta = f'{nombre}_{inicio}_{fin}.fasta'
    
    with open(ruta, 'w') as file:
        file.write(seq)


@st.cache_data
def busca_en_arch(opcion_seleccionada:str, elemento:str, min:int, max:int):
    busqueda = ''
    cromosomas_afectados = set()

    patron = re.compile(rf'^{re.escape(elemento)}\t')
    patron_cromosomas_afectados = re.compile(r'^([^\t]+)')
    patron_pos = re.compile(r'^(?!#)[^\t]+\t([^\t]+)')
    #patron_pos1 = re.compile(r'^(?!#)[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)')
    #patron_pos1 = re.compile(r'^(?!#)[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)')
    lineas_deseadas = set(indices_enfermedades_clndn[opcion_seleccionada])
    with open(RUTA, 'r') as arch:
        for indice, linea in enumerate(arch, 1):
            if indice in lineas_deseadas:
                match_cromosomas = patron_cromosomas_afectados.match(linea)
                # Indica el cromosoma afectado por la enfermedad
                if match_cromosomas:
                    cromosomas_afectados.add(match_cromosomas.group(1))
                # Captura solo el cromosoma indicado
                if patron.match(linea):
                    match = patron_pos.match(linea)
                    if int(match.group(1)) >= min and  int(match.group(1)) <= max:
                        busqueda += linea
                    #if int(match.group(1))>max: break

    return list(cromosomas_afectados), busqueda


@st.cache_data
def busca_crom_en_arch(opcion_seleccionada:str):

    cromosomas_afectados = set()

    patron_cromosomas_afectados = re.compile(r'^([^\t]+)')
    #patron_pos = re.compile(r'^(?!#)[^\t]+\t([^\t]+)')
    #patron_pos1 = re.compile(r'^(?!#)[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)')
    #patron_pos1 = re.compile(r'^(?!#)[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)')
    lineas_deseadas = set(indices_enfermedades_clndn[opcion_seleccionada])
    with open(RUTA, 'r') as arch:
        for indice, linea in enumerate(arch, 1):
            if indice in lineas_deseadas:
                match_cromosomas = patron_cromosomas_afectados.match(linea)
                # Indica el cromosoma afectado por la enfermedad
                if match_cromosomas:
                    cromosomas_afectados.add(match_cromosomas.group(1))


    return list(cromosomas_afectados)

    
def descargar(ruta_archivo):
    with open(ruta_archivo, 'r') as f:
        btn = st.download_button(
            label="Descargar archivo VCF",
            data=f,
            file_name=ruta_archivo,
            mime='text/vcf'
        )

def crea_directori_si_no_existe(dir_path):
    """Create directory if it does not exist."""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def elimina_barra(opcion_seleccionada):
    opcion=''
    for caracter in opcion_seleccionada:
        if caracter in ['/']:
            opcion += '['
        else:
            opcion +=caracter
    return opcion




def obtener_rangos(cromosoma_formateados, email)->dict:
    dic_crom_limites = {}
    for elemento in cromosoma_formateados:
        st.toast(f'Consultando en NCBI el rango del cromosoma {elemento}', icon = '💻')
        min_max = fst.fetch_genomic_dna_fastq(elemento, email)
        st.toast(f'Recibida consulta del NCBI del cromosoma {elemento}', icon='✅')
       
        rango = fst.encontrar_rango_bases(min_max)
        dic_crom_limites[elemento]=rango
    return dic_crom_limites

@st.cache_data
def crea_muestra_dic_rangos(nombre_cromosomas):
    nom_NCBI = fst.convertir_nombre_crom_del_vcf_a_ncbi(nombre_cromosomas)
    st.toast(f'Buscando rangos online...', icon = '🕵️‍♂️')
    dic_rango_aux = obtener_rangos(nom_NCBI,st.session_state.email)


    data = [{'Cromosoma': k, 'Rango Máximo': v} for k, v in dic_rango_aux.items()]
    st.write("## Tabla de Datos")
    st.dataframe(data)
    return dic_rango_aux

def seleccion_rango_busqueda_cadena(longitud:int):
    dynamic_col1, dynamic_col2 = st.columns(2)
    
    with dynamic_col1:
        min = st.number_input(label = f'Inicio', min_value=1, max_value=int(longitud),
                                                             value=1,
                                                             )
    with dynamic_col2:
        max = st.number_input(label = f'Final', min_value=1, max_value=int(longitud),
                                                             value=int(longitud),
                                                             )
    return (min,max)


def seleccion_rango_grupos(grupo, lista,cromosoma):
    dynamic_col1, dynamic_col2= st.columns(2)
  
    lista = [int(elemento) for elemento in lista]
    lista = sorted(lista)
    lista = [str(elemento) for elemento in lista]
    with dynamic_col1:
        min = st.number_input(label = f'Inicio de {grupo}', min_value=1, max_value=int(lista[0]),
                                                             value=int(lista[0])-10000,
                                                             )
    with dynamic_col2:
        max = st.number_input(label = f'Final de {grupo}', min_value=int(lista[-1]), max_value=int(dic_rango_crom[cromosoma]),
                                                             value=int(lista[-1])+10000,
                                                                 )
    return (min,max)

def seleccion_rango_cromosoma(seleccionado):
    dic_cons_seq = {}
    dynamic_col1, dynamic_col2= st.columns(2)
    with dynamic_col1:
        min = st.number_input(label = f'Inicio de {seleccionado}', min_value=1, max_value=int(dic_rango_crom[seleccionado]),
                                                             value=1,
                                                             )
    with dynamic_col2:
        max = st.number_input(label = f'Final de {seleccionado}', min_value=1, max_value=int(dic_rango_crom[seleccionado]),
                                                             value=int(dic_rango_crom[seleccionado]),
                                                             )
        
    
    return (min, max)

def seleccion_rangos(seleccionados):
    dic_cons_seq = {}
    dynamic_col1, dynamic_col2, dynamic_col3, dynamic_col4 = st.columns(4)

    for i, elem in enumerate(seleccionados):
        position = i % 2
        if position == 0:
            with dynamic_col1:
                min = st.number_input(label = f'Inicio de {elem}', min_value=1, max_value=int(dic_rango_crom[elem]),
                                                             value=1
                                                             )
            with dynamic_col2:
                max = st.number_input(label = f'Final de {elem}', min_value=1, max_value=int(dic_rango_crom[elem]),
                                                             value=int(dic_rango_crom[elem]),
                                                             )
            dic_cons_seq[elem] = (min, max)
        else:
            with dynamic_col3:
                min = st.number_input(label = f'Inicio de {elem}', min_value=1, max_value=int(dic_rango_crom[elem]),
                                                             value=1,
                                                             )
            with dynamic_col4:
                max = st.number_input(label = f'Final de {elem}', min_value=1, max_value=int(dic_rango_crom[elem]),
                                                             value=int(dic_rango_crom[elem]),
                                                             )
            dic_cons_seq[elem] = (min, max)
    return dic_cons_seq

def intervalo_elimina_ruido():
    dynamic_col1, dynamic_col2, dynamic_col3, dynamic_col4 = st.columns(4)
    with dynamic_col1:
        min_value_input_fre = st.number_input(f'Frecuencias: Eliminar tramo desde este valor mínimo', 
                    value=-0.000010,
                    step=0.000001,
                    format="%.6f",
                    key='min_input_fre',
                )
    with dynamic_col2:
        max_value_input_fre = st.number_input(
                    f'Frecuencias: Eliminar tramo desde este valor máximo', 
                    value=0.000010,
                    step=0.000001,
                    format="%.6f",
                    key='max_input_fre',
                )
    with dynamic_col3:
        min_value_input_ent = st.number_input(f'Entropías: Eliminar tramo desde este valor mínimo', 
                    value=-0.000040,
                    step=0.000001,
                    format="%.6f",
                    key='min_input_ent',
                )
    with dynamic_col4:
        max_value_input_ent = st.number_input(
                    f'Entropías: Eliminar tramo desde este valor máximo', 
                    value=0.000040,
                    step=0.000001,
                    format="%.6f",
                    key='max_input_ent',
                )
    if min_value_input_fre >= max_value_input_fre:
        st.error(f'El mínimo seleccionado para FRECUENCIA es {min_value_input_fre} y es más grande que el máximo seleccionado {max_value_input_fre}')
        return 0,0,0,0
    if min_value_input_ent >= max_value_input_ent:
        st.error(f'El mínimo seleccionado PARA ENTROPÍA es {min_value_input_fre} y es más grande que el máximo seleccionado {max_value_input_fre}')
        return 0,0,0,0
    return min_value_input_fre,max_value_input_fre,min_value_input_ent,max_value_input_ent

def comprobar_seleccion_rangos(dic_query):
    for clave, valor in dic_query.items():
        if int(valor[0]) > int(valor[1]):
            st.error(f'En el Cromosoma {clave}, el mínimo {valor[0]} es más grande que el valor máximo {valor[1]} indicado')
            return False
    return  True

def obtener_nom_arch_fasta(directorio):
    return [archivo for archivo in os.listdir(directorio) if archivo.endswith(".fasta")]

@st.cache_data
def listar_archivos_fasta(directorio):
    """
    Lista los nombres y las rutas completas de todos los archivos .fasta en un directorio dado,
    devolviendo un diccionario con el nombre del archivo (sin la extensión .fasta) como clave y
    una tupla con el nombre completo del archivo y su ruta como valor.

    Args:
    directorio (str): Ruta al directorio donde buscar archivos .fasta.

    Returns:
    dict: Diccionario donde cada clave es el nombre del archivo sin la extensión `.fasta`
          y el valor es una tupla (nombre del archivo, ruta completa).
    """
    archivos_fasta = {}  # Diccionario para almacenar los nombres de archivos y sus rutas completas
   
    for archivo in os.listdir(directorio):
        if archivo.endswith(".fasta"):  # Verificar si el archivo termina con '.fasta'
            ruta_completa = os.path.join(directorio, archivo)  # Construir la ruta completa del archivo
            nombre_sin_extension = archivo[:-6]  # Remover la extensión '.fasta' del nombre del archivo
            archivos_fasta[nombre_sin_extension] = (archivo, ruta_completa)  # Añadir al diccionario

    return archivos_fasta

# Función para determinar si un elemento es numérico y convertirlo para la ordenación
def clave_ordenacion(elemento):
    if elemento.isdigit():  # Verifica si el elemento es numérico (como string)
        return (0, int(elemento))  # Prioridad 0 para números, convierte a int
    else:
        return (1, elemento)  # Prioridad 1 para letras

def crear_estructura_directorios_y_subcarpetas(base_dic):
    # Definir la ruta base para 'ANALISIS/ENFER/CROMOSOMA'
    if base_dic == {}:return None
    ruta_base = os.path.join(os.getcwd(), 'ANALISIS', 'ENFER', 'CROMOSOMA')
    
    # Crear la estructura de directorios si no existe
    os.makedirs(ruta_base, exist_ok=True)
    
    # Diccionario para guardar las rutas de las nuevas carpetas creadas
    rutas_creadas = {}
    
    # Crear subcarpetas dentro de 'CROMOSOMA' basadas en las claves del diccionario base_dic
    for clave in base_dic:
        # Crear la subcarpeta con el nombre de la clave
        nombre_carpeta = str(clave)
        ruta_carpeta = os.path.join(ruta_base, nombre_carpeta)
        os.makedirs(ruta_carpeta, exist_ok=True)
        
        # Guardar la ruta de la carpeta creada en el diccionario de rutas
        rutas_creadas[clave] = ruta_carpeta
    
    return rutas_creadas



def explorar_directorio_y_extraer_informacion(directorio):
    # Patrón de expresión regular para extraer la información del nombre de archivo
    patron = r"^(NC_\d{6})_(\d+)_(\d+)_(\d+)\.fasta$"
    print(directorio)
    # Lista para guardar la información extraída
    informacion_archivos = []

    # Explorar todos los subdirectorios y archivos
    for raiz, _, archivos in os.walk(directorio):
        for nombre_archivo in archivos:
            # Verificar si el archivo es un archivo .fasta
            if nombre_archivo.endswith('.fasta'):
                # Intentar extraer la información usando la expresión regular
                match = re.match(patron, nombre_archivo)
                if match:
                    # Extraer las partes del nombre del archivo
                    nc_nombre, grupo, pos_inicial, pos_final = match.groups()
                    # Guardar la información en la lista
                    informacion_archivos.append({
                        "ruta_completa": os.path.join(raiz, nombre_archivo),
                        "nc_nombre": nc_nombre,
                        "grupo": grupo,
                        "posicion_inicial": pos_inicial,
                        "posicion_final": pos_final
                    })
    
    return informacion_archivos




if __name__ == '__main__':

    if not os.path.exists('Output'):
        os.makedirs('Output')

    if 'loaded' not in st.session_state:
        st.session_state.page = PAGES.index(dic['1'])
        st.session_state['loaded'] = False

    # Inicializa el estado de la sesión si es necesario
    if 'email_validado' not in st.session_state:
        st.session_state['email_validado'] = False

        
    if not st.session_state.email_validado:
        st.warning("Para continuar, se necesita un email")
        email = st.text_input("Por favor, introduce tu email", key="email_input")
        if st.button("Validar Email"):
            
            if validar_email(email):
                st.session_state.email = email
                st.session_state.email_validado = True
                st.success("El email es válido.")
                st.write("Email guardado:", email)
                print(email)
                if st.button("Iniciar"):
                    run()

            else:
                st.error("El email no es válido. Por favor, introduce un email correcto.")
            
    else:
        run()