"""Maqueta el informe; los algoritmos, distancias, metricas y SVG provienen de C++.
Uso: python scripts/generar_informe.py (requiere reportlab).
"""
from pathlib import Path
import csv
import re
import xml.etree.ElementTree as ET
from html import escape
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.graphics.shapes import Drawing, Rect, Path as DrawPath, String

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output' / 'pdf'
OUT.mkdir(parents=True, exist_ok=True)
W = A4[0] - 88
INK = colors.HexColor('#203747')
TEAL = colors.HexColor('#147487')
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name='TitleLab', fontName='Helvetica-Bold', fontSize=18, leading=22, textColor=INK, spaceAfter=12))
styles.add(ParagraphStyle(name='SectionLab', fontName='Helvetica-Bold', fontSize=13, leading=17, textColor=TEAL, spaceBefore=8, spaceAfter=8))
styles.add(ParagraphStyle(name='TextLab', fontName='Helvetica', fontSize=10, leading=14, spaceAfter=8))
styles.add(ParagraphStyle(name='SmallLab', fontName='Helvetica', fontSize=8.4, leading=11.5, spaceAfter=6))
styles.add(ParagraphStyle(name='CodeLab', fontName='Courier', fontSize=8.2, leading=11, spaceAfter=7))
styles.add(ParagraphStyle(name='CellLab', fontName='Helvetica', fontSize=8, leading=10))
story = []

def p(text, style='TextLab'):
    story.append(Paragraph(text, styles[style]))

def heading(text):
    p(text, 'SectionLab')

def page(title):
    if story:
        story.append(PageBreak())
    p(title, 'TitleLab')

def code(text):
    p(escape(text).replace('\n', '<br/>'), 'CodeLab')

def read_csv(folder, file):
    with (ROOT / 'resultados' / folder / file).open(encoding='utf-8-sig') as f:
        return list(csv.reader(f))

def table(rows, widths=None, size=8):
    cells = [[Paragraph(escape(str(x)).replace('\n', '<br/>'), styles['CellLab']) for x in row] for row in rows]
    t = Table(cells, colWidths=widths, repeatRows=1, hAlign='LEFT')
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e7f0f2')),
        ('TEXTCOLOR', (0, 0), (-1, 0), INK),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 7), ('RIGHTPADDING', (0, 0), (-1, -1), 7),
        ('TOPPADDING', (0, 0), (-1, -1), 6), ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
        ('LINEBELOW', (0, 0), (-1, 0), 0.7, TEAL),
        ('LINEBELOW', (0, 1), (-1, -1), 0.25, colors.HexColor('#cdd7dd')),
    ]))
    story.append(t)
    story.append(Spacer(1, 9))

def matrix(folder, file):
    rows = read_csv(folder, file)
    rows[0][0] = 'Taxón'
    table(rows, [W / len(rows[0])] * len(rows[0]))

def svg_figure(folder, method, max_height=300):
    """Renderiza el subconjunto SVG emitido por el programa (rect, path M/V/H y text).
    No recalcula el arbol ni altera longitudes: conserva coordenadas y etiquetas.
    """
    root = ET.parse(ROOT / 'resultados' / folder / (method + '.svg')).getroot()
    width, height = float(root.attrib['width']), float(root.attrib['height'])
    d = Drawing(width, height)
    for elem in root.iter():
        tag = elem.tag.split('}')[-1]
        a = elem.attrib
        if tag == 'rect':
            d.add(Rect(0, 0, width, height, fillColor=colors.white, strokeColor=None))
        elif tag == 'path':
            tokens = a['d'].split()
            if len(tokens) != 7 or tokens[0] != 'M' or tokens[3] != 'V' or tokens[5] != 'H':
                raise ValueError('SVG inesperado')
            x, y, y2, x2 = map(float, [tokens[1], tokens[2], tokens[4], tokens[6]])
            path = DrawPath(strokeColor=colors.HexColor(a['stroke']), strokeWidth=2, fillColor=None)
            path.moveTo(x, height-y)
            path.lineTo(x, height-y2)
            path.lineTo(x2, height-y2)
            d.add(path)
        elif tag == 'text':
            d.add(String(float(a['x']), height-float(a['y']), elem.text or '',
                         fontName='Helvetica', fontSize=float(a['font-size']),
                         fillColor=colors.HexColor(a.get('fill', '#172b3a'))))
    factor = min(W/width, max_height/height)
    d.scale(factor, factor)
    d.width, d.height = width*factor, height*factor
    story.append(d)
    story.append(Spacer(1, 7))

def metric_table(folders):
    rows = [['Datos', 'Método', 'RMSE', 'MAE', 'Error máx.', 'Ramas < 0']]
    for folder, label in folders:
        for row in read_csv(folder, 'metricas.csv')[1:]:
            rows.append([label, row[0].upper()] + row[1:])
    table(rows, [112, 62, 82, 82, 85, W-423])

page('Laboratorio 03a: Árboles filogenéticos')
p('Biología Molecular Computacional CB309 · 2026-II · 12 de septiembre de 2026', 'SmallLab')
heading('1. Objetivos y alcance')
p('Se implementaron UPGMA y Neighbor Joining (NJ) en C++17 para reconstruir árboles a partir de matrices de distancias. Se procesaron las cinco cadenas de la guía, se aplicó la corrección de Jukes-Cantor y se visualizaron las ramas calculadas. La calidad se evaluó mediante errores de reconstrucción y comprobaciones de ultrametricidad y aditividad. [1, 2]')
p('La prueba con apellidos usa Díaz Neyra, proporcionado por el estudiante, y Túpac Valdivia, correspondiente al profesor identificado en la guía. No se dispuso de los apellidos de otros alumnos. Se añadió un conjunto ficticio, identificado de forma separada, para ejercitar agrupamientos con seis entradas.')
heading('2. Diseño de la solución')
p('El programa sigue este flujo: lectura y validación de datos → cálculo de distancias → construcción del árbol → exportación SVG/Newick → comparación entre distancias originales y distancias del árbol. Todo ese procesamiento se ejecuta en C++; Python se utiliza únicamente para maquetar este PDF a partir de los archivos generados.')
table([
    ['Componente', 'Decisión aplicada'],
    ['Representación', 'Matriz simétrica de double; nodos con hijos, longitudes, altura y número de hojas.'],
    ['UPGMA', 'Fusión del par de menor distancia; actualización ponderada por cantidad de hojas.'],
    ['Neighbor Joining', 'Selección mediante Q, cálculo de ramas y reducción de la matriz.'],
    ['Visualización', 'SVG generado directamente desde C++; Newick para reutilización.'],
    ['Trazabilidad', 'CSV de alineamientos y matrices, pasos de fusión y métricas por ejecución.'],
], [105, W-105])
heading('Repositorio de las implementaciones')
p('<link href="https://github.com/GinoSebastianD/Molecular" color="#147487">https://github.com/GinoSebastianD/Molecular</link><br/>Carpeta de entrega: <b>lab3a-Filogenia</b>.', 'TextLab')
p('Los SVG son cladogramas: su geometría muestra la topología; las etiquetas indican las longitudes numéricas. Las ramas negativas se conservan y se marcan en rojo. En NJ, el punto de dibujo no representa un ancestro inferido.', 'SmallLab')

page('3. Construcción de los algoritmos')
heading('3.1. UPGMA')
p('Cada taxón inicia como un grupo de tamaño 1 y altura 0. Se escoge el par (i, j) con distancia mínima y se crea un nodo u. Las ecuaciones implementadas siguen el procedimiento de las diapositivas 53-59; el promedio debe ponderarse por el número de hojas para conservar el peso de cada taxón. [2]')
code('h(u) = d(i,j) / 2\nL(i,u) = h(u) - h(i)\nL(j,u) = h(u) - h(j)\nd(u,k) = (|i| d(i,k) + |j| d(j,k)) / (|i| + |j|)')
p('Se reemplazan i y j por u y se repite hasta obtener una raíz. Usar (d(i,k) + d(j,k))/2 sin tamaños sería incorrecto cuando los grupos tienen distinto número de hojas. El resultado impone igual distancia de la raíz a todas las hojas.')
heading('3.2. Neighbor Joining')
p('Con n grupos activos, se calculan las sumas de filas r(i). Se minimiza Q(i,j), equivalente al criterio d(i,j) - u(i) - u(j) expuesto en clase, porque Q es ese criterio multiplicado por n - 2. [2, 3]')
code('r(i) = sum_k d(i,k)\nQ(i,j) = (n - 2) d(i,j) - r(i) - r(j)\nL(i,u) = d(i,j)/2 + (r(i)-r(j))/(2(n-2))\nL(j,u) = d(i,j) - L(i,u)\nd(u,k) = (d(i,k) + d(j,k) - d(i,j))/2')
p('Se reduce la matriz después de cada unión. Cuando quedan dos grupos se conectan mediante su distancia restante; para dibujar y serializar se inserta un nodo en el punto medio de esa arista. Esta subdivisión no modifica las distancias entre hojas ni establece una raíz biológica.')
heading('3.3. Decisiones y costo computacional')
p('Los empates conservan el primer par en el orden activo de la matriz, con tolerancia absoluta de 10<super>-12</super>. Los grupos supervivientes conservan su orden y el grupo nuevo se añade al final. No se recortan las longitudes negativas de NJ, pues hacerlo cambiaría los resultados y su error.')
p('Ambos constructores usan O(n³) tiempo y O(n²) memoria numérica: hay O(n) fusiones, cada una con barridos y actualización cuadráticos. El alineamiento global por par requiere O(L<sub>i</sub>L<sub>j</sub>) tiempo y memoria; para n cadenas de longitud máxima L, el conjunto de alineamientos cuesta O(n²L²). Las trazas textuales se conservan para esta práctica pequeña.')

page('4. Preparación de las cinco cadenas')
table([['ID', 'Secuencia original', 'Longitud'], ['S1', 'ATTGCCATT', 9], ['S2', 'ATGGCCATT', 9], ['S3', 'ATCCAATTTT', 10], ['S4', 'ATCTTCTT', 8], ['S5', 'ACTGACC', 7]], [55, 350, W-405])
p('Las cadenas se transcribieron exactamente de la guía. Como sus longitudes difieren, comparar únicamente posiciones originales introduciría una decisión implícita de alineamiento. Se eligió alineamiento global Needleman-Wunsch por pares con coincidencia +2, sustitución -1 y hueco -2. El desempate del retroceso es diagonal, arriba, izquierda. Esta parametrización es una decisión de implementación, no un alineamiento biológico certificado.')
heading('4.1. Alineamientos obtenidos')
alignments = read_csv('adn', 'alineamientos.csv')
rows = [['Par', 'Cadena 1 alineada', 'Cadena 2 alineada', 'Puntaje', 'L', 'm', 'Huecos']]
for row in alignments[1:]:
    rows.append([row[0], row[1], row[2], row[3], row[5], row[6], row[7]])
table(rows, [55, 120, 120, 55, 40, 40, W-430])
p('El símbolo "-" representa un hueco. Para estimar sustituciones, se excluyen las columnas que contienen un hueco en cualquiera de las dos secuencias. Se define L como el número de columnas base-base comparables y m como las sustituciones observadas en ellas. Por tanto, p = m/L. Los huecos y el largo total del alineamiento también se guardan en el CSV, aunque no integran el denominador de p. [4]')
p('Los alineamientos por pares pueden implicar homologías diferentes y producen comparaciones con 6 a 9 sitios. No constituyen un alineamiento múltiple común. Esta limitación es relevante: la matriz resultante puede violar la desigualdad triangular y no representar un árbol con longitudes no negativas.', 'SmallLab')

page('5. Matrices de distancias y corrección')
heading('5.1. Proporción de sustituciones observadas')
matrix('adn', 'distancias_p.csv')
heading('5.2. Corrección de Jukes-Cantor (Ecuación 1)')
code('d_JC = -0.75 ln(1 - (4/3)p),   0 <= p < 0.75')
p('Se aplica la ecuación solicitada a cada par y se completa la matriz simétrica. El programa rechaza p ≥ 0.75 y valores no finitos: la corrección deja de producir una distancia real finita y no se sustituye por un número arbitrario. El modelo considera sustituciones entre nucleótidos, no diferencias entre letras de apellidos. [1, 2, 5]')
matrix('adn', 'distancias.csv')
p('<b>Ejemplo S1-S2:</b> m = 1, L = 9, p = 1/9. La corrección da -0.75 ln(23/27) = <b>0.120257</b> sustituciones por sitio. Para S3-S5: m = 4, L = 7 y p = 4/7, por lo que d = <b>1.076313</b>. Una distancia corregida puede ser mayor que 1 porque estima sustituciones acumuladas por sitio.')
p('Los conteos de columnas comparables L, diferencias m y huecos se presentan junto a cada alineamiento en la sección 4.1. Todos los pares cumplen p &lt; 0.75. Los valores publicados se redondean a seis decimales.', 'SmallLab')

page('6. Resultado UPGMA para ADN')
svg_figure('adn', 'upgma', 305)
p('Figura 1. Árbol UPGMA de las cinco secuencias; etiquetas en sustituciones por sitio. La geometría no está a escala. Fuente: ejecución de C++ con hueco -2.', 'SmallLab')
table([
    ['Paso', 'Grupos fusionados', 'Distancia', 'Altura nueva'],
    [1, 'S1 + S2', '0.120257', '0.060128'],
    [2, 'S3 + S4', '0.304099', '0.152049'],
    [3, 'S5 + (S1,S2)', '0.411980', '0.205990'],
    [4, '(S3,S4) + (S5,(S1,S2))', '0.633424', '0.316712'],
], [42, 259, 104, W-405])
p('La primera unión corresponde al menor valor de la matriz: S1-S2. Después se agrupan S3-S4. S5 se incorpora al grupo de S1 y S2 usando (d(S5,S1) + d(S5,S2))/2 = 0.411980. En la última actualización se pondera por los tamaños de los grupos, no por el número de nodos internos.')
p('Todas las distancias raíz-hoja valen 0.316712. Por ejemplo, para S1 se suman 0.110722 + 0.145861 + 0.060128, con pequeñas diferencias de redondeo. La distancia reconstruida S1-S5 es 0.411980, frente a 0.635473 en la matriz original. El árbol fuerza ultrametricidad y no reproduce exactamente las observaciones.')
code((ROOT/'resultados/adn/upgma.nwk').read_text().strip())

page('7. Resultado Neighbor Joining para ADN')
svg_figure('adn', 'nj', 305)
p('Figura 2. NJ, presentado desde un nodo de dibujo. La rama roja tiene longitud negativa; se conserva su valor algebraico. Fuente: ejecución de C++ con hueco -2.', 'SmallLab')
table([
    ['Paso', 'Grupos fusionados', 'Q mínimo / cierre', 'Ramas nuevas'],
    [1, 'S2 + S5', '-3.543831', '-0.213071; 0.401557'],
    [2, 'S1 + (S2,S5)', '-2.123766', '0.022849; 0.260773'],
    [3, 'S3 + S4', '-1.082360', '0.071139; 0.232960'],
    [4, 'Conexión final entre grupos', 'd = 0.237081', '0.118541; 0.118541'],
], [38, 224, 118, W-380])
p('NJ selecciona primero S2-S5, aunque S1-S2 tiene menor distancia directa, porque utiliza las distancias a los demás grupos. Con n = 5, r(S2) = 1.132702 y r(S5) = 2.976586, el cálculo de L(S2,u) produce -0.213071.')
p('Esta rama no tiene interpretación como longitud evolutiva negativa real. Señala incompatibilidad de las distancias con una reconstrucción de ramas no negativas. No se corrige silenciosamente a cero: el informe y las métricas describen la salida algebraica original. La representación es un árbol sin raíz biológica; su topología distingue la pareja S2-S5 de la pareja S1-S2 favorecida por UPGMA.')
code((ROOT/'resultados/adn/nj.nwk').read_text().strip())

page('8. Calidad de los resultados y ajustes')
p('Para cada par de hojas se suman las longitudes a lo largo de su camino único en el árbol y se comparan con la distancia de entrada. Se usan los 10 pares no ordenados; se excluye la diagonal. RMSE = raíz del promedio de errores al cuadrado; MAE = promedio de errores absolutos. Las métricas usan precisión double antes de redondear los CSV.')
metric_table([('adn', 'ADN, hueco -2'), ('adn_gap1', 'ADN, hueco -1')])
p('En la configuración principal, NJ reduce el RMSE de 0.271091 a 0.048389, aproximadamente 82.2%. Ese mejor ajuste es <b>algebraico</b> e incluye una rama negativa: no demuestra mayor exactitud biológica. No se dispone de un árbol verdadero ni de suficiente información para evaluar esa exactitud.')
heading('8.1. Compatibilidad de la matriz principal')
p('De los 10 triples, 7 incumplen la condición ultramétrica (las dos distancias mayores deberían coincidir). Además, 4 triples incumplen la desigualdad triangular. Un ejemplo es S1-S2-S5: d(S1,S5) = 0.635473 &gt; d(S1,S2) + d(S2,S5) = 0.308743. Esto impide que la matriz sea exactamente la métrica de un árbol con ramas no negativas.')
p('En 4 de los 5 cuartetos tampoco coinciden las dos mayores sumas de la condición de cuatro puntos; el desvío máximo es 0.256998. Estos diagnósticos explican por qué no se espera reconstrucción exacta mediante UPGMA ni NJ. La tolerancia de estos diagnósticos es 10<super>-9</super>.')
heading('8.2. Sensibilidad al costo del hueco')
p('Se repitió el análisis cambiando únicamente el hueco de -2 a -1, manteniendo coincidencia +2 y sustitución -1. Al abaratar los huecos, se excluyen más columnas. S1-S3 pasa de 2/8 sustituciones comparables a 0/7, y S2-S4 pasa de 3/8 a 0/6. Así, secuencias distintas pueden tener p = 0; esto no significa que sean idénticas.')
p('El RMSE de NJ sube a 0.174185 y aparecen tres ramas negativas en la representación exportada. UPGMA baja su error respecto de una <b>matriz diferente</b>; por ello los RMSE entre configuraciones no permiten decidir cuál alineamiento es biológicamente correcto. Se conserva hueco -2 como configuración principal, declarada antes de esta comparación.')
p('Ajustes implementados: alineamiento explícito de cadenas desiguales, exclusión documentada de huecos, validación del dominio de Jukes-Cantor, desempates reproducibles y conservación de ramas negativas. Un estudio biológico requeriría más sitios, un alineamiento múltiple revisado y una evaluación de soporte; no se inventan valores de bootstrap.', 'SmallLab')

page('9. Prueba con los apellidos disponibles')
table([['ID', 'Apellido paterno', 'Apellido materno', 'Procedencia'], ['A1', 'Díaz', 'Neyra', 'Proporcionado por el estudiante'], ['P1', 'Túpac', 'Valdivia', 'Profesor identificado en la guía']], [36, 100, 110, W-246])
p('Se normalizan mayúsculas, tildes y ñ; se eliminan espacios, apóstrofos y guiones. Se calcula Levenshtein con costo unitario de inserción, eliminación y sustitución. Cada distancia se divide entre la longitud del apellido más largo. Para contemplar apellidos cruzados se toma el mínimo de las dos correspondencias posibles:')
code('e(a,b) = Levenshtein(a,b) / max(|a|,|b|)\nd(A,B) = min(e(Ap,Bp)+e(Am,Bm),\n             e(Ap,Bm)+e(Am,Bp)) / 2')
p('Los sufijos p y m indican paterno y materno. Esta medida describe semejanza ortográfica y queda entre 0 y 1; no es una tasa de sustitución de ADN. No se le aplica Jukes-Cantor. Al ignorar el orden, dos personas con el mismo par de apellidos pueden quedar a distancia cero.')
matrix('apellidos_clase', 'distancias.csv')
svg_figure('apellidos_clase', 'upgma', 140)
svg_figure('apellidos_clase', 'nj', 140)
p('Figuras 3 y 4. A1 = Díaz Neyra; P1 = Túpac Valdivia. Ambos métodos reproducen la única distancia, 0.837500, con dos segmentos de 0.418750 y RMSE = 0. Con dos taxones no hay agrupamientos alternativos que comparar ni evidencia de desempeño relativo.', 'SmallLab')
p('<b>¿Hay consistencia con posibles parentescos?</b> No se puede establecer parentesco a partir de este resultado. La distancia indica poca similitud ortográfica entre los apellidos disponibles. No se aportó información genealógica ni la lista de los demás alumnos; la prueba queda limitada a un estudiante y al profesor.')

page('10. Prueba complementaria: datos ficticios')
p('Las siguientes seis entradas son <b>ficticias</b>, construidas para comprobar apellidos compartidos e inversión paterno/materno. No representan alumnos de la clase ni relaciones familiares conocidas.', 'SmallLab')
table([['ID', 'Paterno', 'Materno']] + read_csv('apellidos_demo', 'etiquetas.csv')[1:], [40, 200, W-240])
svg_figure('apellidos_demo', 'upgma', 340)
p('Figura 5. UPGMA sobre los datos ficticios. Conserva las equivalencias D2-D3 y D4-D6 a distancia cero, conforme a la comparación cruzada. D5 se agrupa con Rojas/Silva y D1 con Díaz/Castro. Son asociaciones por escritura; los datos se diseñaron para producir esos contrastes.', 'SmallLab')
page('10.1. NJ en la demostración ficticia')
svg_figure('apellidos_demo', 'nj', 340)
p('Figura 6. Neighbor Joining sobre las mismas seis entradas ficticias. Las etiquetas D1-D6 corresponden a la tabla de la sección 10. La posición del nodo inicial es una convención de dibujo, por lo que las figuras deben compararse por sus conexiones y longitudes, no por la posición visual de la raíz.', 'SmallLab')
metric_table([('apellidos_demo', 'Apellidos ficticios')])
p('El RMSE es 0.035544 para UPGMA y 0.033472 para NJ; ninguna salida tiene ramas negativas. La matriz y los pasos completos están en resultados/apellidos_demo. La coincidencia de apellidos puede sugerir una hipótesis para investigar, pero tampoco en una muestra real probaría parentesco sin información adicional.')
p('Las distancias cero de D2-D3 y D4-D6 son consecuencia exacta del diseño de la medida. Al invertir los dos apellidos, el mínimo cruzado permanece en cero. Esto verifica que el programa considera ambos órdenes y no trata el apellido paterno como la única señal de similitud.')

page('11. Verificación y reproducción')
p('Se compiló con g++ 10.3.0 (TDM-GCC-64), C++17, optimización -O2 y advertencias -Wall -Wextra -Wpedantic. La compilación final terminó sin advertencias. Las <b>19 pruebas automáticas</b> finalizaron correctamente; su registro se entrega en resultados/pruebas.txt.')
table([
    ['Comprobación', 'Evidencia'],
    ['Corrección de distancias', 'p = 0; ejemplo de clase p = 0.25; rechazo de saturación y proporciones negativas.'],
    ['Alineamiento', 'S1-S2 con una sustitución; caso conocido con una deleción.'],
    ['UPGMA', 'Recuperación exacta de una matriz ultramétrica; ejemplo de clase con grupos de tamaños distintos.'],
    ['NJ', 'Recuperación exacta de una matriz aditiva de cinco hojas; cierre con dos taxones; preservación de rama negativa.'],
    ['Robustez', 'Empates deterministas; matrices asimétricas, no cuadradas o negativas rechazadas.'],
    ['Apellidos', 'Normalización; Levenshtein conocido; inversión paterno/materno y simetría.'],
], [112, W-112])
heading('Ejecución desde la carpeta lab3a-Filogenia')
code('powershell -File scripts/reproducir.ps1\n# Linux/macOS con g++:\nsh scripts/reproducir.sh')
p('Los scripts compilan, ejecutan las pruebas y regeneran los cuatro conjuntos: ADN principal, sensibilidad, apellidos disponibles y demostración ficticia. Para cambiar la muestra real, edite data/apellidos_clase.csv o pase otro archivo al script. Los identificadores deben ser únicos; el CSV usa encabezado id,paterno,materno y codificación UTF-8.')
code('g++ -std=c++17 -O2 src/main.cpp -o build/filogenia\n./build/filogenia dna data/secuencias.fasta resultados/adn\n./build/filogenia surnames data/apellidos_clase.csv resultados/apellidos_clase')
p('En Windows se usa build/filogenia.exe. Para regenerar únicamente el informe: instale reportlab y ejecute python scripts/generar_informe.py después de reproducir los resultados. El generador de PDF no implementa UPGMA, NJ ni las distancias.')
heading('Archivos entregados')
p('src/: código C++; tests/: pruebas; data/: entradas; resultados/: matrices, alineamientos, diagnósticos, trazas, árboles SVG y Newick; scripts/: reproducción y maquetación; output/pdf/: informe; README.md: instrucciones. No se incluyen ejecutables ni archivos temporales de compilación.')

page('12. Conclusiones y referencias')
p('1. Se implementaron y ejecutaron UPGMA y NJ en C++17, con visualización propia en SVG. Las pruebas con matrices de referencia verifican la ponderación de UPGMA y la recuperación de distancias aditivas mediante NJ.')
p('2. En las cinco secuencias, UPGMA une primero S1-S2 y NJ une primero S2-S5. NJ presenta menor error respecto de la matriz principal, pero una rama negativa y las violaciones triangulares impiden interpretar la salida como una filogenia de longitudes evolutivas válidas sin revisar las distancias.')
p('3. El alineamiento y el manejo de huecos afectan fuertemente los resultados de estas cadenas cortas. Documentarlos es parte del experimento; aplicar la corrección sin definir sitios comparables ocultaría una decisión esencial.')
p('4. Con Díaz Neyra y Túpac Valdivia ambos métodos ajustan exactamente la única distancia, un caso trivial de dos hojas. La demostración ficticia comprueba el tratamiento de apellidos cruzados, pero no sustituye una comparación con toda la clase ni permite afirmar parentesco.')
heading('Referencias')
p('[1] Túpac Valdivia, Y. J. (2026). <i>Laboratorio 03a: Árboles Filogenéticos</i>. Biología Molecular Computacional CB309, Universidad Católica San Pablo, 9 de septiembre. Documento proporcionado: 03a_UPGMA-NJ_Algorithms_2026-II (1).pdf, pp. 1-2.', 'SmallLab')
p('[2] Túpac Valdivia, Y. J. (2026). <i>Reconstrucción de árboles filogenéticos</i>. Material de clase. Documento proporcionado: CMB_05_Phylogeny_2026-II.pdf. Base de implementación: diapositivas 53-59 (UPGMA), 60-63 (NJ) y 66-67 (corrección de sustituciones).', 'SmallLab')
p('[3] Saitou, N., y Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. <i>Molecular Biology and Evolution</i>, 4(4), 406-425. <link href="https://doi.org/10.1093/oxfordjournals.molbev.a040454" color="#147487">doi:10.1093/oxfordjournals.molbev.a040454</link>.', 'SmallLab')
p('[4] MEGA. <i>Alignment Gaps and Sites with Missing Information</i>. Documentación oficial sobre tratamiento de huecos. <link href="https://www.megasoftware.net/web_help_12/Alignment_Gaps_and_Sites_with_Missing_Information.htm" color="#147487">Consultar documentación</link>. Consulta: 12 de septiembre de 2026.', 'SmallLab')
p('[5] MEGA. <i>Jukes-Cantor Distance Failed</i>. Dominio de aplicación de la corrección. <link href="https://www.megasoftware.net/web_help_12/Jukes-Cantor_Distance_Failed.htm" color="#147487">Consultar documentación</link>. Consulta: 12 de septiembre de 2026.', 'SmallLab')
p('Repositorio: <link href="https://github.com/GinoSebastianD/Molecular/tree/main/lab3a-Filogenia" color="#147487">github.com/GinoSebastianD/Molecular/tree/main/lab3a-Filogenia</link>.', 'SmallLab')

def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor('#ccd7dd'))
    canvas.line(44, 40, A4[0]-44, 40)
    canvas.setFont('Helvetica', 8)
    canvas.setFillColor(INK)
    canvas.drawString(44, 27, 'CB309 · Laboratorio 03a · Implementación en C++')
    canvas.drawRightString(A4[0]-44, 27, str(doc.page))
    canvas.restoreState()

dest = OUT / 'Informe_Laboratorio_03a_Filogenia.pdf'
doc = SimpleDocTemplate(str(dest), pagesize=A4, rightMargin=44, leftMargin=44,
                        topMargin=40, bottomMargin=53, title='Laboratorio 03a: Árboles filogenéticos',
                        author='Díaz Neyra', subject='UPGMA y Neighbor Joining en C++17')
doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(dest)
