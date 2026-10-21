"""Distance conversion and CopraRNA's original tree-weight calculation."""
from pathlib import Path
import re


def transform_distance(path: Path) -> str:
    matrix, organisms = {}, []
    width = None
    for index, line in enumerate(path.read_text().splitlines()[8:],1):
        fields = [v.strip() for v in line.split('\t') if re.search(r'\d',v)]
        if width is None:
            width = len(fields)-1
        organisms.append(re.split(r'\s',fields.pop())[0].lower())
        fields = ['0.00']*(width-len(fields))+fields
        matrix[index] = {j+1:value for j,value in enumerate(fields)}
    # Lexical traversal matches Perl's "sort keys", including >=10 organisms.
    for i in sorted(matrix,key=str):
        for j in sorted(matrix[i],key=str):
            if float(matrix[i][j]) == 0:
                matrix[i][j] = matrix[j][i]
    output = [f'    {width}\n']
    for i in sorted(matrix):
        values = ['0.000' if i==j else matrix[i][j].ljust(5,'0') for j in sorted(matrix[i])]
        output.append(organisms[i-1]+'\t'+''.join(v+'  ' for v in values)+'\n')
    return ''.join(output)


def compute_weights(treefile: Path, neighborfile: Path) -> str:
    tree_length = 0.
    for line in treefile.read_text().splitlines():
        for token in re.split(r'[:(),]',line):
            if re.fullmatch(r'\d+\.\d+',token):
                tree_length += float(token)
    organism_count = 0
    active = False
    rows = []
    for line in neighborfile.read_text().splitlines():
        match = re.search(r'(\d+)\sPopulations',line)
        if match:
            organism_count = int(match[1])
        if re.search(r'From\s+To\s+Length\s+Height',line):
            active = True
        if active and re.search(r'\d',line):
            rows.append(line.replace('-','').strip().split())
    parents = {row[0] for row in rows if row[0].isdigit()}
    children = {row[1] for row in rows if row[1].isdigit()}
    root = next((row[0] for row in reversed(rows) if row[0] in parents-children),'1000000000000000')
    rows.append(['999',root,'0','0'])
    tree, flipped = {}, {}
    for a,b,c,*_ in rows:
        tree.setdefault(a,{})[b] = float(c)
        flipped.setdefault(b,{})[a] = float(c)
    leaves = [child for children in tree.values() for child in children if re.search(r'n[zc]',child)]
    def subtree(node):
        result = 0.
        for child,length in tree.get(node,{}).items():
            result += length
            if not re.search(r'n[zc]',child):
                result += subtree(child)
        return result
    lengths = {row[0]:subtree(row[0]) for row in rows}
    def calculate(root, branch):
        result = 1.
        for nextroot,rootlength in flipped.get(root,{}).items():
            nextsubtree = lengths[root]+rootlength
            if tree_length:
                result = (branch+rootlength/2)/nextsubtree * calculate(nextroot,nextsubtree)
            else:
                result = 1/organism_count
        return result
    weights = {}
    for leaf in leaves:
        for root,branch in flipped[leaf].items():
            for _ in flipped.get(root,{}):
                weights[leaf] = calculate(root,branch)
    return ''.join(f'{leaf};{weights[leaf]:.7f}\n' for leaf in sorted(weights))
