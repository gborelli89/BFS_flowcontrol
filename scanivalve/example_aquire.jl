# Exemplo de aquisição de dados utilizando PyCall

using PyCall

# incluindo caminho da pasta para o pyimport
push!(pyimport("sys")["path"], pwd());
# ou
push!(pyimport("sys")["path"], "/home/Documentos/Doutorado/scanivalve-master");

# Abrindo conexão e configurando
scani = pyimport("scanivalve")
s = scani.Scanivalve(ip="191.30.80.130")
s.config(AVG=1, PERIOD=1000, FPS=1000) # deve dar 16 segundos de medida

# Lendo dados
p,f = s.acquire();
