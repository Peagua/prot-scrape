import pandas as pd
import glob,os
from Bio import SeqIO

input_dir = "blast_filtrado/"
path_proteoma = "proteoma_leish/Leishmania_amazonensis_AnnotatedProteins.fasta"
output_dir = "fastas_hits/"

# Carregando o proteoma
proteoma = SeqIO.to_dict(SeqIO.parse(path_proteoma, "fasta"))

# Função para extrair as sequências, a partir do ID (neste caso ID do Trypdb, presente no proteoma carregado)
def busca_id_proteoma(id, proteoma_usado):
    
    # Busca diretamente o id no proteoma
    if id in proteoma_usado:
        return proteoma_usado[id]

    # Se não encontrar o id, tenta encontrar parte dele (id do trypdb pode aparecer como sujo ou parte do orig)
    for chave in proteoma_usado:
        if id in chave or chave in id:
            return proteoma_usado[chave]

    return None # Caso não encontre de jeito nenhum

top_melhores = int(input("Quantas proteínas quer retornar como melhores hits? "))

arquivos = sorted(glob.glob(f"{input_dir}*.csv"))

for arquivo in arquivos:

    nome = os.path.basename(arquivo).replace("_blast_filtrado.csv", "")
    df = pd.read_csv(arquivo)

    df["score_triagem"]=(
        df["identidade"]*0.5 + df["cobertura_query"]*0.3 + (df["bitscore"]/df["bitscore"].max())*100*0.2
    ) 

    df_ordenado = df.sort_values("score_triagem", ascending=False)

    print(f"\n{nome} - Top {top_melhores}")
    print(df_ordenado[[
        "id_query", "id_limpo_hit", "identidade",
        "cobertura_query", "bitscore", "score_triagem"
    ]].head(top_melhores).to_string(index=False))

    sequencias = []

    for id_hit in df_ordenado["id_limpo_hit"].head(top_melhores):

        seq = busca_id_proteoma(id_hit, proteoma)
        sequencias.append(seq)

    fasta_saida = f"{output_dir}{nome}_melhores_hits.fasta"
    SeqIO.write(sequencias, fasta_saida, "fasta")
    print(f"    {top_melhores} Melhores FASTAS salvos em {fasta_saida}")