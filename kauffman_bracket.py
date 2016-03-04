"""
RelatorCodes is a list of codes of the knot group's relators.
The format of the code is [i,j,k,sign], i is the underpath
 into i-th crossing,j is the underpath from i-th crossing,
 k is the overpath above i-th crossing and sign is the signal
 of i-th crossing.
 
get_StatCrossing takes "RelatorCodes" as an argument and returns
 "Stat_Crossings" which contains the number, it's absolute value is
 the number of a crossing point and it's signals as the crossing.
 So i in Stat_Crossings, i>0 means the i-th crossing is +1,
 i<0 means (-i)-th crossing is -1.
 
toggle_StatCrossings takes 2 arguments, "Stat_Crossing" and "Gauss_Code",
 returns new "Stat_Crossing". This function is needed for the transforming
 to Loo, the orientation of the overpath is changed to the opposite.
 If L(original) is a knot, this reversing doesn't effect the signals
 of crossings, but if L is a link, this reversing effects the crossings
 of the component which has the reversing overpath.
"""
def get_StatCrossings(RelatorCodes):
    Stat_Crossings = []
    for i in RelatorCodes:
        Stat_Crossings.append(i[0]*i[3])
    return(Stat_Crossings)

def toggle_StatCrossings(Stat_Crossings, Gauss_Code):
    newStatCrossings = []
    list_crossings = list(set(map(abs,Gauss_Code)))
    ordr = map(abs, Stat_Crossings)
    for i in Stat_Crossings:
        if abs(i) in list_crossings:
            newStatCrossings.append(-i)
        else:
            newStatCrossings.append(i)
    return(newStatCrossings)

def number2Position(n, stage):
    Position = ''
    bg = Integer(n).digits(2)
    m = len(bg)
    for i in bg:
        Position = str(i)+Position
    for j in range(0, stage - m):
        Position = '0' + Position
    return(Position)
    
"""
Diagram is a list of a GaussCodes and a Stat_Crossings.
Stat_Crossing is a list of integers. Each number is 
a crossing point number and it's sign is a crossing sign.
GaussCodes is a list of GaussCode. 
GaussCode is a list of integers. this list determins 
the oriented knot/link.

Diagram:    [ GaussCodes, Stat_Crossings ]
GaussCodes: [ GaussCode_1,...,GaussCode_m ]
GaussCodes: [ a1, ... a_n ]
         a_i in ZZ, n is a crossing number of the link L.
         a_i<0 means the the path is underpath into the i-th crossing point.
         a_i>0 means the the path is the overpath of i-th crossing point.

To calcukate the Kauffman bracket polynomial of the link L:

    1. Generate the extended gauass code of the link L.
    2. From the extended gauss code, generate a "RelatorCodes" and "GaussCode".
    3. For the Kauffman bracket polynomial, needs only status of crossings,
       which contains a "in-path", a "out-path", a "over-path" and it's signal
       at each crossing points. The information of the path is included
       in "GaussCode", the information of signals of the crossing points is
       included in "Stat_Crossings".
    4. The crossing change from L to Lo and Loo needs the signals 
       of the crossing point. If the signal is 1, then Lo is given 
       by connecting the in-underpath and the out-overpath, 
       the in-overpath and the out-underpath, Loo is given by connecting
       the in-underpath and the in-overpath(orientation reversing),
       the out-overpath(orientation reversing) and the out-underpath.
       If the underpaths and the overpath are in the same component of L, 
       the signals of the crossing points don't change, but those are
       not in the same, the signals of the crossings between the component
       of the underpaths and the ceomponent of the overpath.
       By these crossing change, the given diagram are transformed to 
       list of 2 diagrams Lo, Loo. Lo's factor is A, Loo's factor is A^(-1).
 
        * +1 : [i, j, k, 1]

          L                         Lo                    Loo
          ^                       ^                           |
         j|                      j|                           |k
     k ------->         <=>   ----   ----->         <-------  ------>  
          |                          |k                  l |
         i|                          |                     |
                                                     Orientation reversing 
                                                     on the overpath k.

        * -1 : [i, j, k, -1]

          L                         Lo                    Loo
          ^                      ^                            ^
         j|                     j|                            |j
      <------- k        <=>   ----  ----->          <-------  ------  
          |                         |                      |
         i|                         |k                    i|
                            Orientation reversing 
                            on the overpath k.

5. All crossing changes, there are Diagrams, which GaussCode
   and Stat_Crossings are empty:
       
       Diagrams = [Diagram_1,....Diagram_n]

       Diagram_i : [ [ [] ,..., [] ] , [] ]
       n is 2^(#crossing points=N).
6. Count [] in each entities of Diagrams, #[]-th power of (-A^2 - A^(-2)).
   Diagrams = [Diagram_0, Diagram_1, ... , Diagram_(n-1)]
    Integer      0           1                n-1
    .digits(2) [0,...,0]  [0,...,1], ...., [1,....,1]
    sum()         0          1                N-1
    Coefficient   A^N,    A^N-1*(A^(-1)^1,...,A^N
       
""" 
def crossingChangeA(connectedDiagram):
    [GaussCodes, Stat_Crossings] = connectedDiagram
    GaussCode = GaussCodes[0]
    n = Stat_Crossings[0]
    s = sign(n)
    i = s * n
    y = []
    if len(Stat_Crossings)>1:
        y = Stat_Crossings[1:]
    sa = copy(y)
    sb = copy(y)
    x = []
    p0 = GaussCode.index(-i)
    p1 = GaussCode.index(i)
    if p1>p0:
        "splitting"
        a1 = GaussCode[p1+1:] + GaussCode[0:p0]
        a2 = GaussCode[p0+1:p1]  
        "fusion"
        x = GaussCode[p1-1:p0:-1]
        sb = toggle_StatCrossings(sb, x)
        a0 = [x + a1]
    else:
        "splitting"
        a1 = GaussCode[p1+1:p0]
        a2 = GaussCode[p0+1:]+GaussCode[:p1]
        "fusion"
        x = GaussCode[:p0:-1]
        if p1>0:
            x = GaussCode[p1-1::-1] + x
        sb = toggle_StatCrossings(sb, x)
        a0 = [x + GaussCode[p1+1:p0]]
    if s==1:
        pa = [[a1, a2], sa]
        ma = [a0, sb]
    else:
        pa = [a0, sb]
        ma = [[a1, a2], sa]
    return([pa, ma])

def crossingChangeB(disjointDiagram):
    [GaussCodes, Stat_Crossings] = disjointDiagram
    n = Stat_Crossings[0]
    s = sign(n)
    i = s * n
    y = []
    z = []
    if len(Stat_Crossings)>0:
        y = Stat_Crossings[1:]
    sa = copy(y)
    sb = copy(y)
    if -i in GaussCodes[0] and i in GaussCodes[1]:
        CodeA = GaussCodes[0]
        CodeB = GaussCodes[1]
    else:
        CodeA = GaussCodes[1]
        CodeB = GaussCodes[0]
    p0 = CodeA.index(-i)
    p1 = CodeB.index(i)
    pa =[[CodeB[p1+1::] + CodeB[0:p1] + CodeA[p0+1::] + CodeA[0:p0]], sa]
    x = CodeB[:p1:-1]
    if p1>0:
        x = [CodeB[p1-1::-1] + x]
    sb = toggle_StatCrossings(sb, x)
    pb =[[x + CodeA[p0+1::] + CodeA[0:p0]], sb]
    if s==1:
        z = [pa, pb]
    else:
        z = [pb, pa]
    return(z)

def crossingChange(Diagram):
    [GaussCodes, Stat_Crossings] = Diagram
    i = abs(Stat_Crossings[0])
    y = []
    cpl = []
    newGaussCodes = []
    Diagrams = []
    for x in GaussCodes:
        if not(-i in x or i in x):
            newGaussCodes.append(x)
        else:
            if -i in x and i in x:
                y = crossingChangeA([[x], Stat_Crossings])
            else:
                cpl.append(x)
    if len(cpl)>0:
        y = crossingChangeB([cpl, Stat_Crossings])
    Diagrams = [[y[0][0] + newGaussCodes, y[0][1]],[y[1][0] + newGaussCodes, y[1][1]]]
    return(Diagrams)

def insertDiagramsTable(Rolfsen, Stage, Crossing, Diagrams):
    cursor = sqlite3.connect("/Users/yokotahiroshi/My_KNOT.db")
    m = 0
    for i in Diagrams:
        Position = number2Position(m, Stage)
        m = m + 1 
        GaussCodes = str(i[0])
        Stat_Crossings = str(i[1])
        cursor.execute("insert into diagrams values (?,?,?,?,?,?)",(Rolfsen, int(Stage), int(Crossing), Position, GaussCodes, Stat_Crossings))
    cursor.commit()
    cursor.close()
    
def insertKauffmanBracketTable(Rolfsen, Crossing, Diagram, Polynomial):
    cursor = sqlite3.connect("/Users/yokotahiroshi/My_KNOT.db")
    GaussCode = str(Diagram[0])
    Stat_Crossings = str(Diagram[1])
    cursor.execute("insert into kauffman_bracket values (?,?,?,?,?)",(Rolfsen, int(Crossing), GaussCode, Stat_Crossings, str(Polynomial)))       
    cursor.commit()
    cursor.close()

def calc_Kauffman_Bracket_DB(Rolfsen, Diagram):
    var('A')
    As = []
    Fs = []
    KBTK = []
    Stage = 0
    insertDiagramsTable(Rolfsen, Stage, int(0), [Diagram])
    Stage = 1
    Crossing = abs(Diagram[1][0])
    Diagrams = crossingChange(Diagram)
    insertDiagramsTable(Rolfsen, Stage, Crossing, Diagrams)
    cps = Diagrams[0][1]
    for i in cps:
        Stage = Stage + 1
        tmp = map(crossingChange, Diagrams)
        Diagrams = []
        for j in tmp:
            Diagrams = Diagrams + j
        Crossing = abs(i)
        insertDiagramsTable(Rolfsen, Stage, Crossing, Diagrams)
    for i in Diagrams:
        j = len(i[0]) - 1
        KBTK.append((-A**2-A**(-2))**j)
    n = len(KBTK)
    m = len(Integer(n-1).digits(2))
    Polynomial = 0
    for i in range(0,n):
        As.append(sum(Integer(i).digits(2)))
    for i in As:
        Fs.append(A**(m-i)*A**(-i))
    for i in range(0,n):
        Polynomial = Polynomial + KBTK[i]*Fs[i]
    Polynomial = expand(Polynomial)
    insertKauffmanBracketTable(Rolfsen, Crossing, Diagram, Polynomial)
    return(Polynomial)

def calc_Kauffman_Bracket_polynomial_DB(Rolfsen,Diagram):
    Polynomial = calc_Kauffman_Bracket_DB(Rolfsen, Diagram)
    w = sum(map(sign, Diagram[1]))
    return(expand(-A^(-w)*Polynomial))
