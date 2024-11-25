from Bio import Entrez
#import vcf
import math
import re
from datetime import datetime
import os
#from . import vcf41 as vc
import json
import vcfpy
from collections import Counter
from collections import defaultdict
#from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor


def fetch_genomic_dna_fastq(chromosome:str, start='', end='', email='email@tuemail.com'):
    """
    Obtiene una secuencia de ADN en formato FASTQ desde la base de datos de nucleótidos de NCBI.
    
    Este fragmento de código configura el correo electrónico del usuario para cumplir con las políticas de NCBI,
    realiza una solicitud para obtener la secuencia de ADN del cromosoma especificado entre las posiciones de inicio
    y fin, y devuelve la secuencia en formato FASTQ.
    
    Parámetros:
    - chromosome (str): Identificador del cromosoma o secuencia de referencia de NCBI desde el cual obtener la secuencia de ADN.
    - start (int): Posición de inicio en el cromosoma o secuencia de referencia para la descarga de la secuencia.
    - end (int): Posición final en el cromosoma o secuencia de referencia para la descarga de la secuencia.
    - email (str): Correo electrónico del usuario para uso de la API de Entrez. Es necesario según las políticas de NCBI.
    
    Retorna:
    - str: Secuencia de ADN en formato FASTQ obtenida desde NCBI.
    
    Ejemplo de uso:
    secuencia_fastq = fetch_genomic_dna_fastq(chromosome="NC_000001", start=10000, end=50000, email="tu.email@ejemplo.com")
    """
    
    # Establece el correo electrónico del usuario para las solicitudes de NCBI.
    Entrez.email = email
    
    # Realiza la solicitud a NCBI para obtener la secuencia en formato FASTQ.
    handle = Entrez.efetch(db="nucleotide", rettype="fastq", retmode="text",
                           id=chromosome, seq_start=start, seq_stop=end)
    
    # Lee la secuencia obtenida de la respuesta.
    seq_record = handle.read()
    
    # Cierra el manejador de la conexión.
    handle.close()
    
    # Devuelve la secuencia en formato FASTQ.
    return ''.join(seq_record.splitlines())


def encontrar_rango_bases(secuencia):
    """Busca el rango de bases en la secuencia y retorna una tupla (inicio, fin) o None."""
    patron = re.compile(r"LOCUS\s+\S+\s+(\d+)\s+bp")

    # Buscar el patrón en el texto
    resultado = patron.search(secuencia)
    if resultado:
        return resultado.group(1)
    else:
        return None



def busca_base_inicial_y_base_final(lista_cromosoma:list, email='email@tuemail.com'):
    """Por cada elemento de la lista de un cromosoma devuelve una diccionario cuya clave es el cromosoma y su valor es una tupla"""
    dic ={}
    for elemento in lista_cromosoma:
        aux = fetch_genomic_dna_fastq(str(elemento), email=email)
        rango_bases = encontrar_rango_bases(aux)
        if rango_bases:
            print(f"{elemento} --> Inicio de base: {rango_bases[0]}; Fin de base: {rango_bases[1]}")
            dic[elemento] = rango_bases
        else:
            print(f"No se encontró el patrón buscado para {elemento}.")
            dic[elemento] = None
    return dic




def fetch_genomic_dna_archivo(chromosome:str, start:int, end:int, email='your.email@example.com', filename_prefix=''):
    """
    Descarga secuencias de ADN genómico desde la base de datos de nucleótidos de NCBI y las guarda en un archivo en formato FASTA.
    
    Parámetros:
    - chromosome (str): Identificador del cromosoma o secuencia de referencia de NCBI desde el cual obtener la secuencia de ADN.
    - start (int): Posición de inicio en el cromosoma o secuencia de referencia para la descarga de la secuencia.
    - end (int): Posición final en el cromosoma o secuencia de referencia para la descarga de la secuencia.
    - email (str, opcional): Correo electrónico del usuario para uso de la API de Entrez. Es necesario según las políticas de NCBI. Por defecto, 'your.email@example.com'.
    - filename (str, opcional): Nombre del archivo donde se guardará la secuencia de ADN descargada. Por defecto, 'chromosome_sequence.fasta'.
    
    Funcionamiento:
    1. Configura el correo electrónico del usuario para cumplir con las políticas de NCBI.
    2. Realiza una solicitud a la base de datos de nucleótidos de NCBI para obtener la secuencia de ADN especificada.
    3. Lee la secuencia de ADN obtenida y la guarda en un archivo local especificado.
    
    Ejemplo de uso:
    fetch_genomic_dna_archivo(chromosome="NC_000001", start=10000, end=50000, email="tu.email@ejemplo.com", filename="human_chr1_segment.fasta")
    
    Nota:
    Es importante usar un correo electrónico válido y tener en cuenta las restricciones de uso de la API de Entrez para evitar bloqueos.
    """

    # Genera la fecha y hora actuales
    fecha_hora_actual = datetime.now()
    fecha_hora_str = fecha_hora_actual.strftime("%Y%m%d_%H%M%S")
    
    # Construye el nombre del archivo con la fecha y hora actuales
    if filename_prefix:
        filename = f"{filename_prefix}.fasta"
    else:
        filename = f"{chromosome}_{fecha_hora_str}.fasta"
    
    Entrez.email = email

    # Obtener la secuencia de ADN desde NCBI
    handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=chromosome, seq_start=start, seq_stop=end)
    seq_record = handle.read()
    handle.close()
    
    # Dividir la secuencia en líneas y omitir la primera línea (encabezado)
    lines = seq_record.splitlines()[1:]

    # Guarda la secuencia en un fichero, sin saltos de línea entre las secuencias
    with open(filename, 'w') as file:
        for line in lines:
            file.write(line)  
    print(f"Secuencia guardada en {filename}")


def leer_secuencia_de_fasta_file(nombre_de_archivo):
    """Lee un archivo FASTA y devuelve la secuencia de ADN como una cadena."""
    with open(nombre_de_archivo, 'r') as archivo:
        # Concatenar el resto de las líneas para formar la secuencia
        secuencia = ''.join(linea.strip() for linea in archivo)
        #secuencia = archivo.readline().strip()
    return secuencia

def busca_primera_ocurrencia_en_fasta_file(nombre_de_archivo:str, sqb:str)->tuple:
    # Antiguo nombre: buscar_secuencia_en_fasta_file
    """
    Busca una primera ocurrencia dentro de un archivo FASTA.
    Por ejemplo, busca en un archivo de un cromosoma descargado la posición donde se ubica la secuencia
    a buscar spb
    Parámetros:
    - nombre_de_archivo: El nombre del archivo FASTA desde el cual leer la secuencia.
    - secuencia_a_buscar: La secuencia de caracteres que se busca dentro del archivo FASTA.
    
    Retorna:
    - Una tupla (posición inicial, posición final) de la secuencia buscada dentro del archivo FASTA,
      o None si la secuencia no se encuentra.
    """
    with open(nombre_de_archivo, 'r') as archivo:
        linea = archivo.read()

    posicion_inicial = re.search(re.escape(sqb), linea)

    if posicion_inicial is not None:
        return (posicion_inicial.start() + 1, posicion_inicial.end())
    else:
        return None

def buscar_secuencia_en_fasta_file_todo(nombre_de_archivo:str, sqb:str)->list:
        with open(nombre_de_archivo, 'r') as archivo:
            linea = archivo.read()
        patron = re.escape(sqb)
        posiciones = [(match.start()+1,match.end()) for match in re.finditer(patron, linea)]
        return posiciones

def buscar_secuencia_en_fasta_file_todo_acotado(nombre_de_archivo:str, sqb:str, inicio:int, fin:int)->list:
        with open(nombre_de_archivo, 'r') as archivo:
            linea = archivo.read()
        patron = re.escape(sqb)
        posiciones = []
        for match in re.finditer(patron, linea):
            if inicio <= match.start()+1 and fin >= match.end():
                posiciones.append((match.start()+1,match.end()))
            if fin < match.end(): break
        return posiciones

def busca_primera_ocurrencia_en_fasta(secuencia:str, sqb:str)->tuple:
    """
    Busca una secuencia de caracteres dentro de otra secuencia de ADN proporcionada.
    
    Parámetros:
    - secuencia: La secuencia de ADN completa en la cual buscar.
    - secuencia_a_buscar: La secuencia de caracteres que se busca dentro de la secuencia completa.
    
    Retorna:
    - Una tupla (posición inicial, posición final) de la secuencia buscada dentro de la secuencia completa,
      o None si la secuencia no se encuentra.
    """
    posicion_inicial = re.search(re.escape(sqb), secuencia)

    if posicion_inicial is not None:
        return (posicion_inicial.start() + 1, posicion_inicial.end())
    else:
        return None
    
def buscar_secuencia_en_fasta_todo(secuencia:str, sqb:str)->list:
        patron = re.escape(sqb)
        return [(match.start()+1,match.end()) for match in re.finditer(patron, secuencia)]

def buscar_secuencia_en_fasta_acotado(secuencia:str, sqb:str, inicio:int, fin:int)->list:
        patron = re.escape(sqb)
        posiciones = []
        for match in re.finditer(patron, secuencia):

            if inicio <= match.start()+1 and fin >= match.end():
                posiciones.append((match.start()+1,match.end()))
            if fin < match.end(): break
        return posiciones



def extraer_subsecuencia_de_fasta_file(nombre_de_archivo, inicio, fin):
    """
    Extrae una subsecuencia de un archivo FASTA dado un rango específico, con validaciones adicionales.
    
    Parámetros:
    - nombre_de_archivo: El nombre del archivo FASTA desde el cual leer la secuencia completa.
    - inicio: Posición inicial de la subsecuencia deseada (base 0).
    - fin: Posición final de la subsecuencia deseada (incluida en el resultado).
    
    Retorna:
    - La subsecuencia de caracteres encontrada dentro del rango especificado, o None si los parámetros son inválidos o exceden los límites de la secuencia.
    """
    # Leer la secuencia completa del archivo FASTA
    secuencia_completa = leer_secuencia_de_fasta_file(nombre_de_archivo)
    inicio -=1
    fin -=1
    # Verificar que inicio es menor o igual que fin y que ambos están dentro de los límites de la secuencia
    if inicio <= fin and 0 <= inicio < len(secuencia_completa) and 0 <= fin < len(secuencia_completa):
        # Extraer y retornar la subsecuencia
        return secuencia_completa[inicio:fin+1]  # +1 para incluir el carácter en la posición 'fin'
    else:
        # Retornar None si las condiciones no se cumplen
        return None

def extraer_subsecuencia_de_fasta_online(secuencia_completa, inicio, fin):
    """
    Extrae una subsecuencia de un archivo FASTA dado un rango específico, con validaciones adicionales.
    
    Parámetros:
    - secuencia_completa: secuencia completa de tipo string.
    - inicio: Posición inicial de la subsecuencia deseada (base 0).
    - fin: Posición final de la subsecuencia deseada (incluida en el resultado).
    
    Retorna:
    - La subsecuencia de caracteres encontrada dentro del rango especificado, o None si los parámetros son inválidos o exceden los límites de la secuencia.
    """
    # Leer la secuencia completa del archivo FASTA
    inicio -=1
    fin -=1
    # Verificar que inicio es menor o igual que fin y que ambos están dentro de los límites de la secuencia
    if inicio <= fin and 0 <= inicio < len(secuencia_completa) and 0 <= fin < len(secuencia_completa):
        # Extraer y retornar la subsecuencia
        return secuencia_completa[inicio:fin+1]  # +1 para incluir el carácter en la posición 'fin'
    else:
        # Retornar None si las condiciones no se cumplen
        return None



    
#OBSOLETO    
def buscar_secuencia_en_fasta_online(secuencia, secuencia_a_buscar):
    """
    Busca una secuencia de caracteres dentro de otra secuencia de ADN proporcionada.
    
    Parámetros:
    - secuencia: La secuencia de ADN completa en la cual buscar.
    - secuencia_a_buscar: La secuencia de caracteres que se busca dentro de la secuencia completa.
    
    Retorna:
    - Una tupla (posición inicial, posición final) de la secuencia buscada dentro de la secuencia completa,
      o None si la secuencia no se encuentra.
    """
    # Buscar la secuencia_a_buscar dentro de la secuencia proporcionada
    posicion_inicial = secuencia.find(secuencia_a_buscar)
    
    if posicion_inicial != -1:
        # Calcular la posición final sumando la longitud de la secuencia buscada a la posición inicial
        posicion_final = posicion_inicial + len(secuencia_a_buscar) - 1
        # Ajuste para hacer la posición base-1 en la salida
        return (posicion_inicial + 1, posicion_final + 1)
    else:
        # Devolver None si la secuencia buscada no se encuentra
        return None



def fetch_genomic_dna_online(chromosome:str, start:int, end:int, email='your.email@example.com'):
    """
    Descarga secuencias de ADN genómico desde la base de datos de nucleótidos de NCBI y las guarda en un string.
    
    Parámetros:
    - chromosome (str): Identificador del cromosoma o secuencia de referencia de NCBI desde el cual obtener la secuencia de ADN.
    - start (int): Posición de inicio en el cromosoma o secuencia de referencia para la descarga de la secuencia.
    - end (int): Posición final en el cromosoma o secuencia de referencia para la descarga de la secuencia.
    - email (str, opcional): Correo electrónico del usuario para uso de la API de Entrez. Es necesario según las políticas de NCBI. Por defecto, 'your.email@example.com'.
    
    Funcionamiento:
    1. Configura el correo electrónico del usuario para cumplir con las políticas de NCBI.
    2. Realiza una solicitud a la base de datos de nucleótidos de NCBI para obtener la secuencia de ADN especificada.
    3. Lee la secuencia de ADN obtenida y la guarda en un archivo local especificado.
    
    Ejemplo de uso:
    fetch_genomic_dna_archivo(chromosome="NC_000001", start=10000, end=50000, email="tu.email@ejemplo.com", filename="human_chr1_segment.fasta")
    
    Nota:
    Es importante usar un correo electrónico válido y tener en cuenta las restricciones de uso de la API de Entrez para evitar bloqueos.
    """

    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text",
                           id=chromosome, seq_start=start, seq_stop=end)
    seq_record = handle.read()
    handle.close()
    return ''.join(seq_record.splitlines()[1:])




def remove_first_line_and_newlines(file_path, modified_file_path):
    """
    Elimina la primera línea de un archivo de texto y todos los saltos de línea dentro del archivo.
    
    Parámetros:
    - file_path: Ruta al archivo original.
    - modified_file_path: Ruta al archivo modificado (puede ser el mismo para sobrescribir).
    """
    # Abre el archivo original para leer
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Lee todas las líneas y omite la primera

    # Abre el archivo modificado para escribir
    with open(modified_file_path, 'w') as modified_file:
        for line in lines:
            modified_file.write(line.strip())  # .strip() elimina los espacios en blanco y saltos de línea al principio y al final de cada línea

def generar_diccionario_kmers(secuencia, k):
    """
    Genera un diccionario de k-mers dada una secuencia de ADN y un valor de k.
    
    Parámetros:
    - secuencia (str): La secuencia de ADN de la cual extraer los k-mers.
    - k (int): El tamaño de los k-mers a generar.
    
    Retorna:
    - Un objeto Counter donde las claves son los k-mers únicos encontrados en la secuencia
      y los valores son las veces que cada k-mer se encuentra en la secuencia.
    """
    # Crear diccionario kmers
    dic_fre_kmers = defaultdict(int)
    for i in range(len(secuencia) - k + 1):
        dic_fre_kmers[secuencia[i:i+k]]+=1

    return dic_fre_kmers

def contar_multiplicidades(diccionario_kmers):
    """
    Transforma el diccionario de k-mers en un diccionario de multiplicidades.
    
    Parámetros:
    - diccionario_kmers: El diccionario original de k-mers.
    
    Retorna:
    - Un diccionario donde las claves son las multiplicidades y los valores son
      el número de k-mers que tienen esa multiplicidad.
    """
    multiplicidades = {}
    for kmer, cuenta in diccionario_kmers.items():
        if cuenta in multiplicidades:
            multiplicidades[cuenta] += 1
        else:
            multiplicidades[cuenta] = 1
    return multiplicidades




def contar_caracteres(cadena):
    """Cuenta las ocurrencias de cada carácter en una cadena usando defaultdict para mayor eficiencia."""
    conteo_caracteres = defaultdict(int)  # Inicializa el conteo en 0 automáticamente para cada nueva clave
    for caracter in cadena:
        conteo_caracteres[caracter] += 1
    return dict(conteo_caracteres)


def calcular_probabilidad(dic):
    """Calcula la probabilidad de cada carácter en el diccionario."""
    suma_de_valores = sum(dic.values())
    dic_probab = {}
    for clave, conteo in dic.items():
        dic_probab[clave] = conteo / suma_de_valores
    return dic_probab



def hamming_distance(s1, s2):
    """
    Calcula la distancia de Hamming entre dos cadenas de igual longitud.
    """
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def findMotif(motif, seq, d):
    """
    Busca el motivo de ADN dentro de la secuencia de ADN permitiendo un número específico de mutaciones.
    Devuelve un diccionario con la posición y el motivo encontrado.
    """
    motifs_found = {}
    motif_length = len(motif)
    
    for i in range(len(seq) - motif_length + 1):
        current_substr = seq[i:i+motif_length]
        if hamming_distance(motif, current_substr) <= d:
            motifs_found[i+1] = current_substr
    
    return motifs_found


def Hmax(resultado):
    """Calcula la entropía máxima dada la distribución de probabilidad."""
    tam_dic = len(resultado)
    return sum(-(1/tam_dic) * math.log2(1/tam_dic) * valor for clave, valor in resultado.items())

def Hactual(resultado, prob):
    """Calcula la entropía actual dada la distribución de probabilidad."""
    suma_total = 0
    for clave, probabilidad in prob.items():
        if probabilidad > 0:
            termino = -probabilidad * math.log2(probabilidad)
            suma_total += termino
    return suma_total

def transform_element(element):
    if element.isdigit():  
        return f"NC_0000{element.zfill(2)}"
    elif element == 'X':
        return "NC_000023"
    elif element == 'Y':
        return "NC_000024"
    elif element == 'MT':
        return "NC_012920"
    else:
        return element  
    
def obtener_rango_online(chromosome:str, email='email@tumeail.com'):
    max_1 = fetch_genomic_dna_fastq(chromosome, email)
    return encontrar_rango_bases(max_1)  

def convertir_nombre_crom_del_vcf_a_ncbi(ruta:str)->set:
    conjunt_crom = obtain_all_names_chrom_inVCF(ruta)
    cromosoma_formateados = {transform_element(chrom) for chrom in conjunt_crom}
    return cromosoma_formateados


def leer_secuencia_de_fasta_file_acotada(nombre_de_archivo, inicio=None, fin=None):
    """
    Lee un archivo FASTA y devuelve un segmento de la secuencia de ADN como una cadena.
    
    Parámetros:
    - nombre_de_archivo: El nombre del archivo FASTA desde el cual leer la secuencia.
    - inicio: La posición inicial desde la cual empezar a leer la secuencia (base 1). Si es None, se empieza desde el principio.
    - fin: La posición final hasta la cual leer la secuencia (inclusive). Si es None, se lee hasta el final de la secuencia.
    
    Retorna:
    - Un segmento de la secuencia de ADN contenida en el archivo.
    """
    with open(nombre_de_archivo, 'r') as archivo:      
        # Leer y concatenar el resto de las líneas para formar la secuencia completa
        secuencia_completa = ''.join(linea.strip() for linea in archivo)
        
    # Ajustar los índices para la extracción de la subsecuencia, teniendo en cuenta la indexación base 1
    if inicio is None:
        inicio = 0
    else:
        inicio -= 1  # Ajuste porque el parámetro inicio es base 1, mientras que la indexación en Python es base 0
    
    if fin is None:
        fin = len(secuencia_completa)
    
    # Extraer la subsecuencia especificada
    return secuencia_completa[inicio:fin]




def agrupar_kmers_por_multiplicidad(diccionario_kmers):
    """
    Agrupa k-mers por su multiplicidad de manera eficiente usando defaultdict.
    
    Parámetros:
    - diccionario_kmers: El diccionario original de k-mers.
    
    Retorna:
    - Un diccionario donde las claves son las multiplicidades y los valores son
      listas de k-mers que tienen esa multiplicidad.
    """
    # Utilizar defaultdict para evitar comprobaciones de existencia de clave
    kmers_por_multiplicidad = defaultdict(set)
    for kmer, multiplicidad in diccionario_kmers.items():
        kmers_por_multiplicidad[multiplicidad].add(kmer)
    
    return dict(kmers_por_multiplicidad)





def procesar_kmers(secuencia:str, k:int):
    """
    Procesa una secuencia para obtener k-mers y construir tres diccionarios:
    1. Uno que mapea cada k-mer a su multiplicidad.
    2. Otro que mapea la multiplicidad a la cantidad de k-mers que tienen esa multiplicidad.
    3. Otro que mapea la multiplicidad a los k-mers que tienen esa multiplicidad.
    4. Otro que mapea cada k-mer a su frecuencia relativa.
    
    Parámetros:
    - secuencia: La cadena de caracteres (string) de la secuencia de ADN.
    - k: El tamaño del k-mer.
    
    Retorna:
    - Un diccionario donde las claves son k-mers y los valores son sus multiplicidades.
    - Un diccionario donde las claves son las multiplicidades y los valores son la cantidad de k-mers.
    - Un diccionario donde las claves son las multiplicidades y los valores son listas de k-mers.
    - Un diccionario donde las claves son k-mers y los valores son sus frecuencias relativas.
    """
    # Contar los k-mers en la secuencia
    diccionario_kmers = defaultdict(int)
    for i in range(len(secuencia) - k + 1):
        kmer = secuencia[i:i+k]
        diccionario_kmers[kmer] += 1
    
    total_kmers = len(secuencia) - k + 1  # Total de k-mers en la secuencia

    # Diccionarios de resultados
    multiplicidad_a_cantidad = defaultdict(int)
    multiplicidad_a_kmers = defaultdict(list)
    frecuencias_relativas = {}

    # Llenar todos los diccionarios en un solo paso
    for kmer, multiplicidad in diccionario_kmers.items():
        multiplicidad_a_cantidad[multiplicidad] += 1
        multiplicidad_a_kmers[multiplicidad].append(kmer)
        frecuencias_relativas[kmer] = multiplicidad / total_kmers  # Calcular la frecuencia relativa
    
    return dict(diccionario_kmers), dict(multiplicidad_a_cantidad), dict(multiplicidad_a_kmers), dict(frecuencias_relativas)


def distancia_entre_grupos(rangos:list,distancia:int):
    rangos_inte = [int(numero) for numero in rangos]

    # Lista que contendrá todas las listas de números según las diferencias
    listas_de_numeros = []

    # Lista actual donde se agregan los números
    lista_actual = []
    for i in range(len(rangos_inte)):
    # Si es el primer elemento, solo se agrega a la lista actual
        if i == 0:
            lista_actual.append(rangos_inte[i])
        else:
            # Calcular la diferencia con el número anterior
            diferencia = rangos_inte[i] - rangos_inte[i - 1]
            
            # Si la diferencia es mayor a 10000, crear una nueva lista
            if diferencia > distancia:
                # Agregar la lista actual a la lista de listas y empezar una nueva
                listas_de_numeros.append(lista_actual)
                lista_actual = [rangos_inte[i]]
            else:
                # Si no, seguir agregando a la lista actual
                lista_actual.append(rangos_inte[i])

    # No olvidar agregar la última lista al finalizar el recorrido
    listas_de_numeros.append(lista_actual)

    conversion = [[str(elemento) for elemento in lista] for lista in listas_de_numeros]
    return conversion

def convertir_nombre_crom_del_vcf_a_ncbi(conjunt_crom:set)->set:
    cromosoma_formateados = {transform_element(chrom) for chrom in conjunt_crom}
    return cromosoma_formateados

# Alberto
def obtain_all_names_chrom_inVCF(ruta:str)->set:
    """
    Obtiene un conjunto único de nombres o identificadores desde un archivo VCF
    ubicado en la ruta especificada, basado en una columna específica.

    Args:
        ruta (str): La ruta del archivo VCF del cual se extraerán los nombres o
                    identificadores.


    Returns:
        set: Un conjunto de cadenas únicas que representan los nombres o
             identificadores extraídos de la columna especificada del archivo VCF.

    Nota:
        Este método depende de la ejecución de comandos de shell específicos de Unix.
    """
    pattern = re.compile(r'^(?!#)([^\t]+)')
    with open(ruta,'r') as arch:
        conjunto = {pattern.match(linea).group(1) for linea in arch if pattern.match(linea)}

    return conjunto

def extraer_info_archivo(nombre_archivo):
    """
    Extrae el cromosoma, la posición de inicio y fin del nombre de un archivo FASTA.
    
    Args:
    nombre_archivo (str): El nombre del archivo, ej. 'NC_000003_10000000_11000000.fasta'
    
    Returns:
    dict: Un diccionario con el cromosoma, inicio, y fin extraídos.
    """
    # Elimina la extensión .fasta
    nombre_sin_extension = nombre_archivo.rsplit('.', 1)[0]
    
    # Divide el nombre en partes desde el final para asegurar la correcta extracción de inicio y fin
    partes = nombre_sin_extension.rsplit('_', 2)
    
    if len(partes) < 3:
        raise ValueError("El formato del nombre del archivo no es el esperado.")
    
    cromosoma = partes[0]  # El resto del string es el nombre del cromosoma
    inicio = partes[1]
    fin = partes[2]
    
    return {
        'cromosoma': cromosoma,
        'inicio': int(inicio),  # Convierte a int para manejar como número
        'fin': int(fin)         # Convierte a int para manejar como número
    }
def secuencia_modificada(secuencia_original:str, posiciones:list,dic_acciones:dict)->str:
    resultado = []
    primer_elemento = 0
    for tupla in posiciones:
        resultado.append(secuencia_original[primer_elemento:tupla[0]-1])
        
        resultado.append(dic_acciones.get(tupla))
        primer_elemento = tupla[1]


    resultado.append(secuencia_original[primer_elemento:])
    return ''.join(resultado)

def busca_solapamientos_posiciones(dic11: dict)->dict:
    """
    Esta función analiza un diccionario de posiciones y secuencias, y detecta
    conflictos de solapamientos entre las secuencias en base a sus posiciones.
    
    Args:
    dic11 (dict): Un diccionario cuyas claves son posiciones (como str) y los valores
                  son listas con tuplas que contienen secuencias y subsecuencias.

    Returns:
    dict: Un diccionario de conflictos donde las claves son las tuplas de la posición
          mínima y máxima de los conflictos y los valores son listas de posiciones conflictivas.
    """
    pos = None
    dic_conflictos = {}
    conflictos = []
    clave_anterior = None

    for clave, valor in dic11.items():
        clave_int = int(clave)
        longitud = len(valor[0][0])
        fin_actual = clave_int + (longitud - 1)  # Calcula la posición final actual

        if pos is None:  # Primera iteración
            pos = fin_actual
            clave_anterior = clave_int if longitud > 1 else None
            continue

        if clave_int > pos:  # No hay solapamiento
            if conflictos:
                dic_conflictos[(min(conflictos), max(conflictos))] = conflictos
                conflictos = []
            pos = fin_actual
            clave_anterior = clave_int if longitud > 1 else None
        else:  # Hay solapamiento
            if clave_anterior:
                conflictos.append(clave_anterior)
                clave_anterior = None

            conflictos.append(clave_int)
            pos = max(pos, fin_actual)

    # Registra los últimos conflictos si existen
    if conflictos:
        dic_conflictos[(min(conflictos), max(conflictos))] = conflictos

    return dic_conflictos

def modificar_fasta(archivo_entrada, archivo_salida, dic_acciones):
    with open(archivo_entrada, 'r') as entrada, open(archivo_salida, 'w') as salida:
        # Lee la secuencia completa, que es la única línea en el archivo
        secuencia = entrada.read().strip()
        
        # Ordena las acciones de sustitución para procesarlas secuencialmente
        acciones_ordenadas = dic_acciones.items() #sorted(dic_acciones.items())
        
        # Variable para acumular la secuencia resultante
        resultado = []
        primer_elemento = 0
        
        # Procesa cada sustitución en la secuencia
        for clave, valor in acciones_ordenadas:
            inicio, fin = clave
            inicio -= 1  # Ajustar a índice basado en 0
            
            # Agrega parte de la secuencia antes de la sustitución
            resultado.append(secuencia[primer_elemento:inicio])
            resultado.append(valor[1])
            # Actualiza el punto de inicio para la siguiente parte de la secuencia
            primer_elemento = fin
        
        # Agrega el resto de la secuencia después de la última sustitución
        resultado.append(secuencia[primer_elemento:])
        
        # Escribe la secuencia modificada en el archivo de salida
        salida.write(''.join(resultado))

def modificar_fasta_solo(archivo_entrada, dic_acciones):
    with open(archivo_entrada, 'r') as entrada:
        # Lee la secuencia completa, que es la única línea en el archivo
        secuencia = entrada.read().strip()
        
        # Ordena las acciones de sustitución para procesarlas secuencialmente
        acciones_ordenadas = dic_acciones.items() #sorted(dic_acciones.items())
        
        # Variable para acumular la secuencia resultante
        resultado = []
        primer_elemento = 0
        
        # Procesa cada sustitución en la secuencia
        for clave, valor in acciones_ordenadas:
            inicio, fin = clave
            inicio -= 1  # Ajustar a índice basado en 0
            
            # Agrega parte de la secuencia antes de la sustitución
            resultado.append(secuencia[primer_elemento:inicio])
            resultado.append(valor[1])
            # Actualiza el punto de inicio para la siguiente parte de la secuencia
            primer_elemento = fin
        
        # Agrega el resto de la secuencia después de la última sustitución
        resultado.append(secuencia[primer_elemento:])
        
        # Escribe la secuencia modificada en el archivo de salida
        return ''.join(resultado) 
    

def comparar_diccionarios(dic1:dict, dic2:dict)->dict:
    # Nuevo diccionario para almacenar los resultados
    resultado = {}

    # Iterar sobre cada k-mer en el primer diccionario
    for kmer in dic1:
        if kmer in dic2:
            # Calcula la diferencia en valor absoluto si el k-mer está en ambos diccionarios
            resultado[kmer] = abs(dic1[kmer] - dic2[kmer])
        else:
            # Si el k-mer no está en el segundo diccionario, toma el valor en negativo
            resultado[kmer] = -dic1[kmer] 

    # Verifica los k-mers en el segundo diccionario que no fueron procesados
    for kmer in dic2:
        if kmer not in dic1:
            # Agrega el valor en negativo de estos k-mers al resultado
            resultado[kmer] = -dic2[kmer] 

    return resultado