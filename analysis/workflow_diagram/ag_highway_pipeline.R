# ########################################################## ####
# AG introductions in KLEBSIELLA                             ####
# Author:    Cara Conradsen                                  ####
# ########################################################## ####



# AG pipeline ----------------------------------------------------------
DiagrammeR::grViz("digraph{
graph [layout = dot, fontname = Arial, rankdir = TB]

node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25, fontname = Arial]
kpne_assemblies[label =< Hybrid assemblies from E. Feil <br/> <i> K.pneumoniae, n = 1,695</i>>]
pirate[label = <Get pangenome<br/><font face='Courier New'>PIRATE</font>>]
get_ags[label = <Get AG loci from <font face='Courier New'>kpne_ag </font>project>]

subgraph cluster_0 {
graph[shape = rectangle]
bgcolor = grey90
color = grey90

label = \"nBLAST databases\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
gtdb[label = <Genome Taxonomy Database<br/><font face='Courier New'>GTDB 226</font>>]
uhgg[label = <Unified Human Gastrointestinal Genome<br/><font face='Courier New'>UHGG v2.0.2</font>>]
}

kpne_assemblies -> pirate -> get_ags -> {gtdb uhgg}
}")



