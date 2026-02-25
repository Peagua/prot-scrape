# Pipeline para encontrar proteínas em Leishmania

Projeto para encontrar proteínas potencialmente interativas com certas moléculas (é pro meu TCC).
Muitas estruturas nesse projeto seguem padrões para utilizar na minha pesquisa.
Não necessáriamente é reprodutível para todos tipos de dados.

## Workflow
Inicialmente, deve-se procurar o smiles da molécula desejada no PubChem, e buscar nas abas de "Identity" e "Similarity" pela opção de "Proteins" em "Linked Data". Essa opção mostra proteínas alvo dos estudos envolvendo a molécula identidade (se existir) ou similiares ao smiles forncecido.
Em seguida, deve-se baixar o .csv das proteínas, contendo as informações de "protid", "Taxonomy" e "Protein_Accession". Após baixar os arquivos, deve-se preparar eles para o script get_sequence_ref.
    - Dentro de "proteinas_triagem", criar um diretório com o nome identificador do composto (Ex.: RMS)
    - Mover os arquivos baixados para a pasta criada
    - Os arquivos devem seguir um identificador, se a origem foi de moléculas identidade ou similar (Ex.: RMS_I.csv)
A partir dai, é só executar o script, alterando a lista de compostos com os nomes dos identificadores de cada csv de cada molécula testada.