#import vcf
import vcfpy
import re
import subprocess
import matplotlib.pyplot as plt
import json
#import pprint
import time
import orjson

# Patrón de expresión regular para extraer el texto deseado
patron = r"##INFO=<ID=([^,]+),"
pattern_pos = re.compile(r'^(?!#)[^\t]*\t([^\t]+)')
#pp =  pprint.PrettyPrinter(indent=4)
elementos_escalpar = ['.','^', '$', '*', '[', ']', '+', '?', '{', '}', '\\','(',')']


def escapar_caracteres(cadena):
    """ 
    Lee el contenido de la cadena y si encuentra caracteres especiales añade delante la \ 
    para que la expresión regular se tome en cuenta ese caracter  como literal y no como elemento de la expresión regular.
    """
    cadena_escapada = ""
    for caracter in cadena:
        if caracter in elementos_escalpar:
            cadena_escapada += '\\' + caracter
        else:
            cadena_escapada += caracter
    return cadena_escapada


# Alberto
def encuentra_param_vcf4(nombre_archivo):
    """
    Lee el contenido de un archivo y encuentra todas las coincidencias de un patrón definido.

    La función abre y lee un archivo cuyo nombre se pasa como argumento. Elimina los espacios
    en blanco y los saltos de línea de cada línea leída, y luego busca coincidencias de un 
    patrón predefinido en la secuencia de texto resultante. Devuelve una lista de todas las 
    coincidencias encontradas.

    EXTRAE todos los SUBCAMPOS del campo FIJO de INFO. 

    Args:
        nombre_archivo (str): El nombre o ruta del archivo que se va a leer.

    Returns:
        list: Una lista de todas las coincidencias del patrón en el archivo. Si no se encuentra
              ninguna coincidencia, devuelve una lista vacía.

    Notas:
        - La variable 'patron' debe ser definida fuera de la función y debe contener una
          expresión regular válida.
        - Importa el módulo  're' para utilizar esta función correctamente.
    """
    with open(nombre_archivo, 'r') as archivo:
        # Combina todas las líneas del archivo en una sola cadena, eliminando espacios
        # en blanco y saltos de línea.
        secuencia = ''.join(linea.strip() for linea in archivo)

    # Busca todas las coincidencias del patrón en la secuencia de texto.
    resultados = re.findall(patron, secuencia)

    return resultados



# Alberto
def lector_vcf(vcf_path):

    """
    Crea un objeto lector para leer datos de un archivo VCF.

    Esta función abre un archivo VCF ubicado en la ruta especificada y crea un objeto lector
    utilizando la biblioteca 'vcf'. Esto permite acceder a los datos contenidos en el archivo
    VCF de manera estructurada, facilitando la lectura de los mismos.

    Args:
        vcf_path (str): La ruta al archivo VCF que se desea leer. Puede ser una ruta absoluta
                        o relativa.

    Returns:
        vcf.Reader: Un objeto lector asociado al archivo VCF especificado. Este objeto puede
                    ser utilizado para iterar sobre los registros del archivo VCF y acceder a
                    la información contenida en cada uno de ellos.

    Notas:
        - Es NECESARIO asegurar que el archivo VCF exista en la ruta especificada antes
          de llamar a esta función para evitar errores de archivo no encontrado.
        - Esta función asume que la biblioteca 'vcf' está instalada y disponible en el entorno
          de ejecución. La biblioteca 'vcf' debe proporcionar, al menos, un método o clase
          llamado 'Reader' capaz de leer archivos VCF.
        - El archivo se abre en modo de lectura ('r'). Asegúrate de cerrar el objeto lector
          cuando haya terminado de usarlo para liberar recursos del sistema.
    """
    return  vcfpy.Reader(open(vcf_path, 'r'))

# Alberto
def numero_datos_vcf(nombre_archivo):
    """
    Cuenta el número de registros de datos en un archivo VCF, excluyendo las líneas de metadatos.

    Esta función abre un archivo VCF y lee su contenido línea por línea e ignora cualquier línea
    que empiece con el carácter '#' pues son líneas de metadatos en el formato VCF. Sólo
    considera para el conteo aquellas líneas que representan registros de datos, es decir, 
    aquellas que no son parte de los metadatos del archivo. Esto permite obtener el número de 
    registros de variantes genéticas presentes en el archivo.

    Args:
        nombre_archivo (str): La ruta al archivo VCF que se va a analizar. Este puede ser
                              un camino absoluto o relativo.

    Returns:
        int: El número total de registros de datos encontrados en el archivo VCF, excluyendo
             las líneas de metadatos.

    Ejemplo:
        Si tienes un archivo VCF con 100 registros de datos y 10 líneas de metadatos, esta 
        función devolverá 100.
    
    Notas:
        - Es importante asegurarse de que el archivo exista en la ruta especificada antes de
          intentar abrirlo para evitar errores.
        - La función automáticamente cierra el archivo después de leer su contenido, gracias
          al uso del gestor de contexto `with`.
    """
    with open(nombre_archivo, 'r') as archivo:
        #ignora las líneas de metadatos
        lineas_datos = [linea for linea in archivo if not linea.startswith('#')]
    return len(lineas_datos)


# Alberto
def buscar_id_en_vcf(ruta:str, numeroID)->dict:
    """
    Busca una ID específica en un archivo VCF y extrae información relacionada.

    El archivo VCF contiene solo estos campos:
        CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    Esta función utiliza el comando `grep` de UNIX para buscar una ID específica en un archivo VCF,
    basándose en una expresión regular construida alrededor de la ID. Si la ID se encuentra, extrae y
    procesa información relacionada a la variante, incluyendo cromosoma, posición, ID, referencia, 
    alteración, calidad, filtro e información adicional. La información adicional (INFO) se procesa 
    adicionalmente para convertirla en un diccionario de pares clave=valor.

    Args:
        ruta (str): La ruta al archivo VCF donde se buscará la ID.
        numeroID (str): La ID específica a buscar dentro del archivo VCF.

    Returns:
        dict: Un diccionario con la información extraída para la ID encontrada, con las claves siendo
              los campos del archivo VCF y los valores siendo los datos correspondientes a esos campos.
              Si se encuentra información INFO adicional, esta se devuelve como un diccionario anidado.
              Retorna None si la ID no se encuentra en el archivo.

    Notas:
        - Esta función utiliza `subprocess.run` para ejecutar `grep` y buscar la ID, por lo que es 
          necesario que `grep` esté disponible en el sistema donde se ejecute la función.
        - La función decodifica la salida de `grep` a UTF-8, asumiendo que el archivo VCF está codificado
          en este formato.
    """
    # Construye la expresión regular para buscar la ID específica en el archivo.
    expresion_regular = f'^([^\t]*\t){{2}}{numeroID}(\t|;)' 
    
    # Ejecuta `grep` para buscar la expresión regular en el archivo.
    resultado = subprocess.run(['grep', '-E', expresion_regular, ruta], stdout=subprocess.PIPE)

    # Decodifica y procesa la salida de `grep`.
    lineas_encontradas = resultado.stdout.decode('utf-8').strip().split('\n')

    if lineas_encontradas[0]:
        print(f'ID {numeroID} encontrada en el archivo VCF:')
        registro_vcf={}
        for linea in lineas_encontradas:

            columnas = linea.split('\t')
            registro_vcf[columnas[2]] = {
            'CHROM': columnas[0],
            'POS': columnas[1],
            'ID': columnas[2],
            'REF': columnas[3],
            'ALT': columnas[4],
            'QUAL': columnas[5],
            'FILTER': columnas[6],
            'INFO': columnas[7]
            }

            # Para la información INFO, que contiene varios pares clave=valor separados por ";", la podemos procesar adicionalmente
            info_items = registro_vcf[columnas[2]]['INFO'].split(';')
            info_dict = {item.split('=')[0]: item.split('=')[1] if '=' in item else item for item in info_items}

            # Agrega la información procesada al diccionario
            registro_vcf[columnas[2]]['INFO'] = info_dict
        return registro_vcf
    else:
        print(f'La ID {numeroID} no se encontró en el archivo.')
        return None

# Alberto
def buscar_pos_en_vcf(ruta:str, numeroPOS)->dict:
    """
    Busca entradas específicas por su número de POS en un archivo VCF y extrae información relevante de ellas.
    El archivo VCF contiene solo estos campos:
        CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

    La función utiliza el comando 'grep' de Unix para buscar en el archivo VCF especificado por la ruta
    las líneas que coincidan con el 'numeroPOS' proporcionado. Después, procesa estas líneas para extraer
    y organizar los datos de interés en un diccionario.

    Args:
        ruta (str): La ruta del archivo VCF en el cual buscar.
        numeroPOS (str): El identificador específico a buscar dentro del archivo VCF.

    Returns:
        dict: Un diccionario que contiene los datos extraídos de las entradas encontradas para el 'numeroID'
              especificado. Cada entrada del diccionario representa una línea del archivo VCF donde el 'ID'
              coincide con el 'numeroID' proporcionado. Si no se encuentra el 'numeroID', devuelve None.

    Notas:
        - Esta función depende del comando 'grep' de Unix, sólo funcionará correctamente en
          sistemas compatibles con Unix.
        - Asume que el archivo VCF está estructurado en columnas separadas por tabuladores (TSV), siguiendo
          el formato estándar VCF.
    """

    # Define una expresión regular para buscar el numeroID específico en el archivo.
    expresion_regular = f'^([^\t]*\t){{1}}{numeroPOS}(\t|;)' 

    # Ejecuta el comando 'grep' para buscar la expresión regular en el archivo.
    resultado = subprocess.run(['grep', '-E', expresion_regular, ruta], stdout=subprocess.PIPE)

    # Decodifica el resultado a texto y lo divide por líneas.
    lineas_encontradas = resultado.stdout.decode('utf-8').strip().split('\n')
    
    if lineas_encontradas[0]:
        print(f'POS {numeroPOS} encontrada en el archivo VCF:')
        registro_vcf={}

        for linea in lineas_encontradas:
            # Divide la línea en columnas basado en tabuladores.
            columnas = linea.split('\t')
            
            # Almacena la información relevante en un diccionario.
            registro_vcf[columnas[2]]= {
            'CHROM': columnas[0],
            'POS': columnas[1],
            'ID': columnas[2],
            'REF': columnas[3],
            'ALT': columnas[4],
            'QUAL': columnas[5],
            'FILTER': columnas[6],
            'INFO': columnas[7]
            }

     
            info_items = registro_vcf[columnas[2]]['INFO'].split(';')
            info_dict = {item.split('=')[0]: item.split('=')[1] if '=' in item else item for item in info_items}

       
            registro_vcf[columnas[2]]['INFO'] = info_dict

        return registro_vcf
    else:
        print(f'La POS {numeroPOS} no se encontró en el archivo.')
        return None
    

# Alberto
def buscar_campo_en_vcf(ruta:str, cadena:str, elemento:str)->dict:
    """
    Busca un elemento específico que contenga un SUBCAMPO del CAMPO INFO dentro de un archivo VCF
    utilizando expresiones regulares y el comando 'grep' de UNIX porque es más eficiente.
    El archivo VCF contiene solo estos campos:
        CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

    SUBCAMPO son todos aquellos  elementos que contiene el CAMPO INFO. Por ejemplo:
    'ALLELEID', 'CLNDISDB', 'CLNDN', 'GENEINFO', 'CLNREVSTAT'...

    Esta función ejecuta un comando 'grep' para buscar en un archivo VCF (Formato de
    Archivo de Variante) específico por un elemento que contenga una cadena dada.
    Construye un diccionario que mapea los IDs de variantes a sus respectivos detalles
    obtenidos de las líneas del archivo VCF que coinciden con el criterio de búsqueda.
    NOTA: IMPORTANTE QUE TODOS LOS ID DEL VCF SEAN ÚNICOS PUES SE UTILIZA COMO CLAVE EN LOS DICCIONARIOS


    Args:
        ruta (str): La ruta al archivo VCF en el que se realizará la búsqueda.
        cadena (str): EL SUBCAMPO de INFO se buscará dentro del valor del 'elemento'
                      especificado en las entradas INFO del archivo VCF.
        elemento (str): El nombre del elemento dentro de la entrada INFO en el archivo VCF
                        que debe contener la cadena de texto especificada.

    Returns:
        dict: Un diccionario que contiene los detalles de las variantes que coinciden con
              la búsqueda. Cada clave del diccionario es el ID de una variante, y su valor
              es otro diccionario con los detalles de esa variante (CHROM, POS, ID, REF,
              ALT, QUAL, FILTER, INFO). El campo INFO se devuelve como un diccionario de
              elementos INFO subdivididos. Si no se encuentran coincidencias, se devuelve None.

    Nota:
        Esta función depende del comando 'grep' de UNIX, por lo que su ejecución está limitada
        a sistemas que dispongan de este comando. Instala 'grep' y que
        el archivo VCF esté correctamente formateado y accesible en la ruta especificada.
    """


    # Eliminar esta línea, al final. Es otra versión
    # expresion_regular = f'{elemento}=([^;]*{cadena}[^;]*)'
    cadena = escapar_caracteres(cadena)
    # Construye la expresión regular para buscar la cadena dentro del valor del elemento especificado.
    expresion_regular = f'{elemento}=(?:[^;]*\|)?({cadena})(?:\||;|$)'

    # Ejecuta el comando 'grep' para buscar la expresión regular en el archivo. 
    resultado = subprocess.run(['grep', '-E', expresion_regular, ruta], stdout=subprocess.PIPE)
  
    # Decodifica la salida de 'grep' y divide el resultado en líneas.
    lineas_encontradas = resultado.stdout.decode('utf-8').strip().split('\n')

    # Inicializa el diccionario que contendrá los registros VCF encontrados.
    registro_vcf ={}

    # Procesa cada línea encontrada que coincida con el criterio de búsqueda.
    if lineas_encontradas[0]:
        print(f'{elemento} {cadena} encontrado en el archivo VCF:')

        for linea in lineas_encontradas:
            # Separa los campos de la línea VCF por tabuladores.
            columnas = linea.split('\t')

            # Construye un diccionario con los detalles de la variante.
            registro_vcf[columnas[2]] = {
            'CHROM': columnas[0],
            'POS': columnas[1],
            'ID': columnas[2],
            'REF': columnas[3],
            'ALT': columnas[4],
            'QUAL': columnas[5],
            'FILTER': columnas[6],
            'INFO': columnas[7]
            }

         
            info_items = registro_vcf[columnas[2]]['INFO'].split(';')
            info_dict = {item.split('=')[0]: item.split('=')[1] if '=' in item else item for item in info_items}

            registro_vcf[columnas[2]]['INFO'] = info_dict
            

        return registro_vcf
    else:
        print(f'El {elemento} {cadena} no se encontró en el archivo.')
        return None
    

# Alberto
def dibuja_grafica_datos(dic, total_datos):
    """
    Genera y muestra gráficas de barras para cada clave en un diccionario, comparando el número
    de datos con y sin un determinado valor (representado por la clave).

    Para cada clave en el diccionario dado, esta función calcula la diferencia entre el total de
    datos y el valor asociado a la clave, generando una comparación visual entre el número de datos
    'Con [clave]' y 'Sin [clave]'. Cada par de categorías se representa en su propia gráfica de
    barras, donde las barras también muestran el porcentaje del total de datos que representan.

    Args:
        dic (dict): Un diccionario donde cada clave es una categoría y el valor asociado es el
                    número de datos que caen en esa categoría.
        total_datos (int): El número total de datos considerados, usado para calcular los
                           datos sin la categoría especificada y los porcentajes.

    Esta función utiliza matplotlib para crear y mostrar las gráficas. Cada gráfica incluye dos
    barras: una para los datos 'Con [clave]' y otra para los datos 'Sin [clave]'. Se añade un
    título y etiquetas a los ejes para mejorar la comprensión de la gráfica. Además, se calcula
    y muestra el porcentaje que cada barra representa respecto al total de datos, posicionando
    el texto del porcentaje encima de cada barra para facilitar su lectura.

    Notas:
        - Es necesario tener instalado y configurado matplotlib para utilizar esta función.
        - La función muestra las gráficas inmediatamente después de crearlas con plt.show().
          Esto puede interrumpir la ejecución del código hasta que se cierre cada ventana de gráfica.
    """

    
    for clave, valor in dic.items():
        datos_sin_clndn = total_datos - valor

        # Nombres de las categorías
        categorias = [f'Con {clave}', f'Sin {clave}']

        # Valores correspondientes a cada categoría
        valores = [valor, datos_sin_clndn]

        # Crear gráfica de barras
        fig, ax = plt.subplots()
        bars = ax.bar(categorias, valores, width=0.5, color=['blue', 'orange'], zorder=3)
        for bar in bars:
            height = bar.get_height()
            porcentaje = f'{(height / total_datos) * 100:.1f}%'
            ax.text(bar.get_x() + bar.get_width() / 2., height, porcentaje, ha='center', va='bottom')
        # Añadir título y etiquetas a los ejes
        ax.set_title('Distribución de datos con y sin CLNDN')
        ax.set_ylabel('Número de líneas')
        ax.set_xlabel('Categoría')

        # Mostrar la gráfica
        plt.show()


# Alberto
def extraer_cabeceras_vcf(nombre_archivo:str)->str:
    """
    Extrae todas las líneas de cabecera de un archivo VCF y las guarda en un string.

    Esta función abre y lee un archivo VCF dado, extrayendo todas las líneas que comienzan
    con el carácter '#', lo cual indica que son cabeceras o comentarios en el formato VCF.
    Cada línea de cabecera se agrega a un string, incluyendo un salto de línea al final
    de cada una para mantener la estructura original del archivo VCF cuando se visualice
    el string o se guarde en un nuevo archivo.

    Args:
        nombre_archivo (str): La ruta y el nombre del archivo VCF a procesar.

    Returns:
        str: Un string que contiene todas las líneas de cabecera del archivo VCF,
        cada una separada por un salto de línea.
    """
    cabeceras = ""
    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            if linea.startswith('#'):
                # Agregamos cada cabecera al string, incluyendo un salto de línea
                cabeceras += linea.strip() + '\n'
    return cabeceras


# Alberto
def buscar_campo_en_vcf_string(ruta:str, cadena:str, elemento:str)->str:
    """
    Busca un elemento específico que contenga un SUBCAMPO del CAMPO INFO dentro de un archivo VCF
    utilizando expresiones regulares y el comando 'grep' de UNIX porque es más eficiente.
    El archivo VCF contiene solo estos campos:
        CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

    SUBCAMPO son todos aquellos  elementos que contiene el CAMPO INFO. Por ejemplo:
    'ALLELEID', 'CLNDISDB', 'CLNDN', 'GENEINFO', 'CLNREVSTAT'...

    Esta función ejecuta un comando 'grep' para buscar en un archivo VCF (Formato de
    Archivo de Variante) específico por un elemento que contenga una cadena dada.
    Construye un diccionario que mapea los IDs de variantes a sus respectivos detalles
    obtenidos de las líneas del archivo VCF que coinciden con el criterio de búsqueda.
    NOTA: IMPORTANTE QUE TODOS LOS ID DEL VCF SEAN ÚNICOS PUES SE UTILIZA COMO CLAVE EN LOS DICCIONARIOS


    Args:
        ruta (str): La ruta al archivo VCF en el que se realizará la búsqueda.
        cadena (str): EL SUBCAMPO de INFO se buscará dentro del valor del 'elemento'
                      especificado en las entradas INFO del archivo VCF.
        elemento (str): El nombre del elemento dentro de la entrada INFO en el archivo VCF
                        que debe contener la cadena de texto especificada.

    Returns:
        str: Un string que contiene los detalles de las variantes que coinciden con
              la búsqueda. Cada clave del diccionario es el ID de una variante, y su valor
              es otro diccionario con los detalles de esa variante (CHROM, POS, ID, REF,
              ALT, QUAL, FILTER, INFO). Si no se encuentran coincidencias, se devuelve None.

    Nota:
        Esta función depende del comando 'grep' de UNIX, por lo que su ejecución está limitada
        a sistemas que dispongan de este comando. Instala 'grep' y que
        el archivo VCF esté correctamente formateado y accesible en la ruta especificada.
    """
    j = time.time()
    cadena = escapar_caracteres(cadena)
    # Construye la expresión regular para buscar la cadena dentro del valor del elemento especificado.
    expresion_regular = f'{elemento}=(?:[^;]*\|)?({cadena})(?:\||;|$)'

    # Ejecuta el comando 'grep' para buscar la expresión regular en el archivo. 
    resultado = subprocess.run(['grep', '-E', expresion_regular, ruta], stdout=subprocess.PIPE)
    
    lineas_encontradas = resultado.stdout.decode('utf-8').strip().split('\n')
    ele = ''
    if lineas_encontradas[0]:
        for linea in lineas_encontradas:
            ele  += linea + '\n'
        h = time.time()
        print(f'tiempo tardado {h-j}')
        return ele
    else:
        print(f'El {elemento} {cadena} no se encontró en el archivo.')
        return None


# Alberto
def guardar_string_como_vcf(contenido:str, nombre_archivo_destino:str)->None:
    """
    Guarda una cadena de texto en un archivo con extensión .vcf.

    Esta función toma una cadena de texto (`contenido`) y un nombre de archivo
    destino (`nombre_archivo_destino`). Asegura que el nombre del archivo tenga
    la extensión .vcf. Si el nombre proporcionado no termina en .vcf, la función
    automáticamente añade esta extensión al nombre del archivo. Posteriormente,
    abre el archivo en modo de escritura y guarda el contenido de la cadena de texto
    en el archivo. Finalmente, imprime un mensaje de éxito indicando que el archivo
    ha sido guardado satisfactoriamente.

    Args:
        contenido (str): La cadena de texto que se desea guardar en el archivo.
        nombre_archivo_destino (str): El nombre del archivo donde se guardará el contenido.
                                      Si el nombre del archivo no termina en .vcf, la función
                                      automáticamente añade esta extensión al nombre.

    Returns:
        None: Esta función no devuelve ningún valor. Imprime un mensaje de éxito al terminar.
    """
    # Asegura de que el nombre del archivo tiene la extensión .vcf
    if not nombre_archivo_destino.endswith('.vcf'):
        nombre_archivo_destino += '.vcf'
    
    # Abre el archivo en modo de escritura y guardar el contenido
    with open(nombre_archivo_destino, 'w') as archivo:
        archivo.write(contenido)
    
    print(f"El archivo '{nombre_archivo_destino}' ha sido guardado satisfactoriamente.")

# Alberto
# def guardar_dicc_json(ruta:str, diccionario:dict)->None:
#     with open(ruta, 'w', encoding='utf-8') as archivo:
#         json.dump(diccionario, archivo, ensure_ascii=False, indent=4)
#     print(f"Diccionario guardado en {ruta}")

# Alberto
def obtain_all_ID_inVCF(ruta:str)->set:
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
    pattern = re.compile(r'^(?!#)[^\t]*\t[^\t]*\t([^\t]+)')
    with open(ruta,'r') as arch:
        conjunto = {pattern.match(linea).group(1) for linea in arch if pattern.match(linea)}

    return conjunto

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

# Alberto
def busca_todos_valores_en_un_subcampo_de_info(ruta:str, subcampo:str)->set:
    """
    Busca todos los valores únicos para un subcampo específico dentro del campo INFO de un archivo VCF.

    Utiliza una expresión regular para identificar líneas que contienen el subcampo especificado,
    extrae esos valores, y retorna un conjunto de valores únicos. Si el valor del subcampo contiene
    múltiples elementos separados por '|', estos se dividen y se tratan como valores separados.

    Parámetros:
    - ruta (str): La ruta al archivo VCF en el que se realizará la búsqueda.
    - subcampo (str): El nombre del subcampo dentro del campo INFO a buscar.

    Retorna:
    - Un conjunto (set) de valores únicos encontrados para el subcampo especificado. Si el subcampo
      no se encuentra en ninguna parte del archivo, retorna None y muestra un mensaje indicativo.

    Ejemplo de uso:
    conjunto_valores = busca_todos_valores_en_un_subcampo_de_info('/path/al/archivo.vcf', 'CLNDN')

    Notas:
    - La función asume un formato VCF donde el campo INFO es la octava columna.
    - Se requiere que 'grep' esté disponible en el entorno donde se ejecuta la función.
    """
    # Construyo la expresión regular para capturar el valor del subcampo especificado
    expresion_regular = f'{subcampo}=([^;]*)'

    # Ejecuta 'grep' para buscar la expresión regular en el archivo
    resultado = subprocess.run(['grep', '-E', expresion_regular, ruta], stdout=subprocess.PIPE)
    lineas_encontradas = resultado.stdout.decode('utf-8').strip().split('\n')

    conjunto = set()

    # Procesa cada línea encontrada para extraer el valor del subcampo
    if lineas_encontradas[0]:  # Verifica que se haya encontrado al menos una línea
        for linea in lineas_encontradas:
            columnas = linea.split('\t')

            # Crea un diccionario a partir del campo INFO
            info_dict = {item.split('=')[0]: item.split('=')[1] if '=' in item else item for item in columnas[7].split(';')}
            
            elem = info_dict.get(subcampo)
            if elem:  
                if '|' in elem:  
                    subelementos = elem.split('|')
                    conjunto.update(subelementos)  
                else:
                    conjunto.add(elem) 
        return conjunto

    else:
        print(f'El {subcampo} no se encontró en el archivo {ruta}.')
        return None
    
def leer_dic_json(ruta_archivo):
    with open(ruta_archivo, 'r') as archivo:
        diccionario_leido = json.load(archivo)
    return diccionario_leido



def buscar_id_en_vcf_str(ruta:str, numeroID)->str:
    """
    Busca una ID específica en un archivo VCF y extrae información relacionada.

    El archivo VCF contiene solo estos campos:
        CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    Esta función utiliza el comando `grep` de UNIX para buscar una ID específica en un archivo VCF,
    basándose en una expresión regular construida alrededor de la ID. Si la ID se encuentra, extrae y
    procesa información relacionada a la variante, incluyendo cromosoma, posición, ID, referencia, 
    alteración, calidad, filtro e información adicional. La información adicional (INFO) se procesa 
    adicionalmente para convertirla en un diccionario de pares clave=valor.

    Args:
        ruta (str): La ruta al archivo VCF donde se buscará la ID.
        numeroID (str): La ID específica a buscar dentro del archivo VCF.

    Returns:
        dict: Un diccionario con la información extraída para la ID encontrada, con las claves siendo
              los campos del archivo VCF y los valores siendo los datos correspondientes a esos campos.
              Si se encuentra información INFO adicional, esta se devuelve como un diccionario anidado.
              Retorna None si la ID no se encuentra en el archivo.

    Notas:
        - Esta función utiliza `subprocess.run` para ejecutar `grep` y buscar la ID, por lo que es 
          necesario que `grep` esté disponible en el sistema donde se ejecute la función.
        - La función decodifica la salida de `grep` a UTF-8, asumiendo que el archivo VCF está codificado
          en este formato.
    """
    # Construye la expresión regular para buscar la ID específica en el archivo.
    expresion_regular = f'^([^\t]*\t){{2}}{numeroID}(\t|;)' 
    
    # Ejecuta `grep` para buscar la expresión regular en el archivo.
    resultado = subprocess.run(['grep', '-E', expresion_regular, ruta], stdout=subprocess.PIPE)

    # Decodifica y procesa la salida de `grep`.
    lineas_encontradas = resultado.stdout.decode('utf-8').strip().split('\n')

    if lineas_encontradas[0]:
        return '\n'.join(lineas_encontradas)
    else:
        print(f'La ID {numeroID} no se encontró en el archivo.')
        return None
    


def buscar_pos_en_vcf_str(ruta:str, numeroPOS)->str:
    """
    Busca entradas específicas por su número de POS en un archivo VCF y extrae información relevante de ellas.
    El archivo VCF contiene solo estos campos:
        CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

    La función utiliza el comando 'grep' de Unix para buscar en el archivo VCF especificado por la ruta
    las líneas que coincidan con el 'numeroPOS' proporcionado. Después, procesa estas líneas para extraer
    y organizar los datos de interés en un diccionario.

    Args:
        ruta (str): La ruta del archivo VCF en el cual buscar.
        numeroPOS (str): El identificador específico a buscar dentro del archivo VCF.

    Returns:
        dict: Un diccionario que contiene los datos extraídos de las entradas encontradas para el 'numeroID'
              especificado. Cada entrada del diccionario representa una línea del archivo VCF donde el 'ID'
              coincide con el 'numeroID' proporcionado. Si no se encuentra el 'numeroID', devuelve None.

    Notas:
        - Esta función depende del comando 'grep' de Unix, sólo funcionará correctamente en
          sistemas compatibles con Unix.
        - Asume que el archivo VCF está estructurado en columnas separadas por tabuladores (TSV), siguiendo
          el formato estándar VCF.
    """

    # Define una expresión regular para buscar el numeroID específico en el archivo.

    expresion_regular = f'^([^\t]*\t){{1}}{numeroPOS}(\t|;)' 

    # Ejecuta el comando 'grep' para buscar la expresión regular en el archivo.
    resultado = subprocess.run(['grep', '-E', expresion_regular, ruta], stdout=subprocess.PIPE)

    # Decodifica el resultado a texto y lo divide por líneas.
    lineas_encontradas = resultado.stdout.decode('utf-8').strip().split('\n')
    
    if lineas_encontradas[0]:
        return '\n'.join(lineas_encontradas)
    else:
        print(f'La ID {numeroPOS} no se encontró en el archivo.')
        return None

def crea_dic_elementos_encontrados(elementos:str)->dict:
    """
        Procesa una cadena de texto que representa elementos de un archivo VCF (Variant Call Format)
        y crea un diccionario con información detallada de cada variante genética encontrada.

        Args:
            elementos (str): Una cadena de texto que contiene líneas del archivo VCF. Cada línea
            representa una variante y sus campos están separados por tabulaciones.

        Returns:
            dict: Un diccionario donde cada clave es el ID de la variante y el valor es otro
            diccionario con detalles de la variante como cromosoma, posición, referencias,
            alternativas, calidad, filtro e información adicional.

        Estructura del diccionario retornado:
        {
            variant_id: {
                'CHROM': str,   # Cromosoma
                'POS': str,     # Posición
                'ID': str,      # ID de la variante
                'REF': str,     # Referencia
                'ALT': str,     # Alternativa
                'QUAL': str,    # Calidad
                'FILTER': str,  # Filtro
                'INFO': dict    # Información adicional procesada en formato diccionario
            },
            ...
        }
    """
    registro_vcf={}
    for linea in elementos:

        columnas = linea.split('\t')
        registro_vcf[columnas[2]] = {
        'CHROM': columnas[0],
        'POS': columnas[1],
        'ID': columnas[2],
        'REF': columnas[3],
        'ALT': columnas[4],
        'QUAL': columnas[5],
        'FILTER': columnas[6],
        'INFO': columnas[7]
        }

        # Para la información INFO, que contiene varios pares clave=valor separados por ";", la podemos procesar adicionalmente
        info_items = registro_vcf[columnas[2]]['INFO'].split(';')
        info_dict = {item.split('=')[0]: item.split('=')[1] if '=' in item else item for item in info_items}

        # Agrega la información procesada al diccionario
        registro_vcf[columnas[2]]['INFO'] = info_dict
    return registro_vcf

def transform_element(element):
    """ 
        Transforma un elemento del nombre del cromosoma que por ejemplo se haya encontrado
        en el VCF y lo cambia a nombre establecido en NCBI
    
    """
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

def guardar_dicc_json(ruta:str,dic:dict):
    with open(ruta, 'wb') as f:
        f.write(orjson.dumps(dic))
    print(f"Diccionario guardado en {ruta}")