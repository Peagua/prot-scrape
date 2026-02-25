import pandas as pd
import requests, time
from tqdm import tqdm
from Bio import Entrez
import shutil, glob
import os


# Input para entrar com o e-mail utilizado para validar o API do NCBI
Entrez.email = str(input("Entre com o email para o Enterz: "))


# Função para padronizar a tabela de proteínas de entrada
def fix_tabel(data):
    # Remove o id de acesso do pubchem CASO NECESSÁRIO
    #data.drop('Protein_Accession', axis=1, inplace=True)

    # Renomeia as colunas
    data = data.rename(columns= {'protid': 'protein', 'Taxonomy': 'organism',
                                 'Protein_Accession': 'pubchem_id'}) 

    # Remove vírgula CASO NECESSÁRIO
    #data["protein"]=data["protein"].str.replace(",", "", regex=False)

    # Removendo o conteúdo do último parenteses (nome do organismo), utilizando expressões regulares
    # r'' para abrir a regex
    # \s* captura o espaço antes do parênteses
    # \( pega o parenteses de abertura
    # [^)]* pega todo conteúdo que não seja ")"
    # *\) pega o fechamento do parenteses
    # \s* pega o espaço após o parenteses
    # $ garante que essa seleção ocorra somente na última ocorrência de parênteses 
    data["protein"]=data["protein"].str.replace(r'\s*\([^)]*\)\s*$', '', regex=True)

    # Na coluna de organismo, pega apenas as 2 primeiras palavras
    data["organism"]=data["organism"].str.split().str[:2].str.join(' ')
    return data

# Função para busca primária, com API do UniProt
def busca_uniprot(nome_proteina, organismo):
    url = "https://rest.uniprot.org/uniprotkb/search" # URL do API do Uniprot

    # Parâmetros de busca no API
    parametros = {
        "query": f'protein_name:"{nome_proteina}" AND organism_name:"{organismo}"',
        "format": "fasta",
        "size": 1
    }

    # Retornando os resultados encontrados
    try:
        # Chamando o API
        resposta = requests.get(url, params=parametros, timeout=10)
        
        # Caso tenha erro do tipo HTTP, se houver
        resposta.raise_for_status()

        # Removendo espaços em branco na busca
        conteudo = resposta.text.strip()
        
        # Todo FASTA começa com >, logo, se a resposta começar com >, retorna o FASTA
        if conteudo.startswith(">"):
            return conteudo
        else:
            return None
    
    # Evidenciando erros
    except requests.exceptions.RequestException as e:
        print(f'\nErro UniProt [{nome_proteina}] | [{organismo}]: {e}')
        return None
    
    # Espera 0.4 segundos, para não dar overload no API (pode aumentar por segurança)
    finally:
        time.sleep(0.4)

# Função de busca usando API do NCBI
def busca_ncbi(nome_proteina, organismo, id_pubchem):
    try:
        if not id_pubchem:
            # Função para encontrar os IDs no NCBI (caso não exista)
            handle_busca = Entrez.esearch(
                db="protein",
                term=f'"{nome_proteina}"[Protein Name] AND "{organismo}"[Organism]',
                retmax=1
            )
            resultado = Entrez.read(handle_busca) #Enterz abre a função e executa
            handle_busca.close() # Fecha a função

            # Pega apenas os ids encontrados, e coloca em uma variável
            ids_encontrados = resultado["IdList"]

            # Caso não encontre nenhum ID, retorna None
            if not ids_encontrados:
                return None
            
            id_usado = ids_encontrados[0] # Apenas o primeiro ID encotrado

        # No caso de ter ID do PubChem
        else:
            id_usado = id_pubchem
        
        # Agora, a partir dos IDs encontrados, realizar a busca pelo Fasta
        handle_fasta = Entrez.efetch(
            db = "protein",
            id=id_usado,
            rettype = "fasta",
            retmode="text"
        )
        fasta = handle_fasta.read().strip() # Limpa o fasta encontrado
        handle_fasta.close() # Como antes, tem que fechar a função

        time.sleep(0.4) # Respeita o tempo do API

        # Caso for um fasta válido (começa com ">"), o retorna
        return fasta if fasta.startswith(">") else None
    
    # Em caso de erro complexo, envolvendo a entrada e/ou API
    except Exception as e:
        if id_pubchem:
            print(f'\nErro NCBI [{nome_proteina} | {organismo} | {id_pubchem}]: {e}')
        else:
            print(f'\nErro NCBI [{nome_proteina} | {organismo}]: {e}')
        time.sleep(1)
        return None

# Função para integrar os mecanísmos de busca, Loopar as procuras, 
# escrever o .fasta com as sequências e relatório das buscas
def recuperar_sequencias(dados, arquivo_saida, comp, busca_id=False):
    # Listas controle, para serem preenchidas com as sequências e relatório
    sequencias_encontradas = []
    relatorio = []

    # Loop usando o tqdm (para mostrar barra de progresso), para executar as buscas em cada proteína
    for _,linha in tqdm (dados.iterrows(), total=len(dados), desc="Buscando Proteínas"):
        
        # variáveis para armazenar nome da proteína e organismo
        nome = linha["protein"]
        org = linha["organism"]

        # variáveis de controle
        fasta = None # para a sequência
        fonte = None # para o relatório

        # Caso deseja realizar a busca pelo ID do Pubchem
        if busca_id == True:
            id_pc = linha["pubchem_id"]
            fasta = busca_ncbi(nome_proteina=nome, organismo=org, id_pubchem=id_pc)
            if fasta:
                fonte = "NCBI"

        if busca_id == False:
            # Primeiro tenta no UniProt
            fasta = busca_uniprot(nome_proteina=nome, organismo=org)
            if fasta:
                fonte = "UniProt"
            
            # Em seguida tenta no NCBI
            if not fasta:
                fasta = busca_ncbi(nome, org, id_pubchem=None)
                if fasta:
                    fonte = "NCBI"
            
        # Caso encontre um fasta, coloca na lista de sequências, e adiciona ao relatório
        if fasta:
            sequencias_encontradas.append(fasta)
            if busca_id == True:
                relatorio.append({
                    "Proteína":  nome,
                    "Organismo": org,
                    "PubChemID": id_pc,
                    "Status":    "encontrada",
                    "Fonte":     fonte
                    })
            
            else:
                relatorio.append({
                    "Proteína":  nome,
                    "Organismo": org,
                    "Status":    "encontrada",
                    "Fonte":     fonte
                    })
        # Se não encontrar, coloca apenas no relatório
        else:
            if busca_id == True:
                relatorio.append({
                    "Proteína":  nome,
                    "Organismo": org,
                    "PubChemID": id_pc,
                    "Status":    "não encontrada",
                    "Fonte":     None
                    })
            
            else:
                relatorio.append({
                    "Proteína":  nome,
                    "Organismo": org,
                    "Status":    "não encontrada",
                    "Fonte":     None
                    })
    
    # Após todas as buscas, escreve o .fasta
    with open(arquivo_saida, "w") as f:
        f.write("\n\n".join(sequencias_encontradas))
    
    # Para criar o dataframe do relatório
    df_relatorio = pd.DataFrame(relatorio)
    df_relatorio.to_csv(f"{comp}_relatorio_busca.csv", index=False)

    total = len(dados)
    encontradas = df_relatorio[df_relatorio["Status"] == "encontrada"].shape[0]
    
    # Printando um apanhado geral no terminal
    print(f'\n Concluído: {encontradas}/{total} proteínas recuperadas')
    if busca_id == False:
        print(f'\n-> UniProt: {df_relatorio[df_relatorio['Fonte']=='UniProt'].shape[0]}')
    
    print(f'\n-> NCBI: {df_relatorio[df_relatorio['Fonte']=='NCBI'].shape[0]}')
    print(f'\n-> Não Encontradas: {total - encontradas}')

    # A função retorna o relatório em dataframe, caso queira printar no terminal
    return df_relatorio


# Input verificador para saber qual tipo de busca realizar
verificador = True
while verificador == True:
    tem_id = str(input("Possui PubChemID nos dados?[S]/[N] "))
    if tem_id == "S":
        usar_id = True
        verificador = False
    elif tem_id == "N":
        usar_id = False
        verificador = False
    else:
        print('Resposta inválida. Tente novamente.')


# Próprio para minha pesquisa, contendo o endereço inicial de cada arquivo utilizado
lista_comps = ['TJL_S','DDS_I','DDS_S','RMS_I','RMS_S']

for i in range(len(lista_comps)):
    composto = lista_comps[i]

    print(f'\nExecutando buscas para o arquivo {composto}.csv')

    # Variável para aramazenar apenas o nome do composto
    nome_composto=composto.replace("_"," ").split()[0]

    # Abrindo a lista com o panda
    tabela = pd.read_csv(f"proteinas_triagem/{nome_composto}/{composto}.csv")

    # Arrumando o df
    tabela = fix_tabel(tabela)

    # Variável para armazenar o nome do arquivo fasta com as sequências
    arquivo = str(f"{composto}_proteinas.fasta")

    # Rodando a função principal
    relat_final=recuperar_sequencias(dados=tabela, arquivo_saida=arquivo, comp=composto, busca_id=usar_id)

    print(relat_final)

# Criando pastas e realocando os fastas e relatórios para as determinadas pastas
os.makedirs('relatorios_busca', exist_ok=True)
for arquivo in glob.glob('*relatorio_busca.csv'):
    shutil.move(arquivo, 'relatorios_busca/')

os.makedirs('fastas_ref', exist_ok=True)
for arquivo in glob.glob('*proteinas.fasta'):
    shutil.move(arquivo, 'fastas_ref/')
