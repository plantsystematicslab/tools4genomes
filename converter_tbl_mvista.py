#!/usr/bin/env python3
import sys
import re

def converter():
    # Verifica se os argumentos de entrada e saída foram passados
    if len(sys.argv) < 3:
        print("Uso: python3 converter.py <arquivo_entrada.tbl> <arquivo_saida.txt>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Erro: O arquivo '{input_file}' não foi encontrado.")
        sys.exit(1)

    output = []
    vistos = set()
    coord_re = re.compile(r'^(\d+)\s+(\d+)\s+(\w+)')

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if not line or line.startswith('>Feature'):
            i += 1
            continue
            
        match = coord_re.match(line)
        if match:
            s_raw, e_raw = int(match.group(1)), int(match.group(2))
            
            # Identifica fita e ordena coordenadas para o mVista
            strand = ">" if s_raw < e_raw else "<"
            start, end = min(s_raw, e_raw), max(s_raw, e_raw)
            
            # Busca o nome do gene/produto nas linhas abaixo
            gene_name = "unknown"
            j = i + 1
            while j < len(lines) and (not lines[j].strip() or lines[j].startswith('\t')):
                attr = lines[j].strip()
                if attr.startswith('gene') or attr.startswith('product'):
                    parts = attr.split()
                    if len(parts) > 1:
                        gene_name = parts[-1]
                j += 1
            
            # Evita duplicar a mesma região (Ex: gene e CDS com mesmas coordenadas)
            key = f"{start}-{end}"
            if key not in vistos:
                output.append(f"{strand} {start} {end} {gene_name}")
                output.append(f"{start} {end} exon\n")
                vistos.add(key)
            
            i = j - 1
        i += 1

    # Grava o resultado no arquivo de saída especificado
    with open(output_file, 'w', encoding='utf-8') as f_out:
        f_out.write("\n".join(output))
    
    print(f"Conversão concluída: '{input_file}' -> '{output_file}'")

if __name__ == "__main__":
    converter()
