# ########################################################## ####
# AG introductions in KLEBSIELLA                             ####
# Author:    Cara Conradsen                                  ####
# ########################################################## ####



# AG pipeline ----------------------------------------------------------
DiagrammeR::grViz("digraph{
graph [layout = dot, fontname = Arial, rankdir = TB]

node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25, fontname = Arial]
pirate[label = <Get pangenome<br/><font face='Courier New'>PIRATE</font>>]
get_ags[label = <Get accessory genes from <font face='Courier New'>kpne_ag </font>project>]
bakta_prot[label =  <annotate proteins <br/><font face='Courier New'> bakta_proteins 1.8.1</font>>]
bakta_anno[label =  <annotate proteins <br/><font face='Courier New'> bakta 1.8.1</font>>]
subset_panx[label =  <subset to <i>n = 485</i>, determine accessory genes>]

subgraph cluster_0 {
graph[shape = rectangle]
bgcolor = grey90
color = grey90

label = \"Pangenome data\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
kpne_assemblies[label =< Long read assemblies from E. Feil <br/> <i> K.pneumoniae, n = 485 </i>>]
panX[label =<panX gene alignments <br/> <i>E.coli </i>&amp; <i>S.aureus</i>>]
}


subgraph cluster_1 {
graph[shape = rectangle]
bgcolor = grey90
color = grey90

label = \"nBLAST databases\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
gtdb[label = <Genome Taxonomy Database<br/><font face='Courier New'>GTDB 226</font>>]
uhgg[label = <Unified Human Gastrointestinal Genome<br/><font face='Courier New'>UHGG v2.0.2</font>>]
}

kpne_assemblies -> bakta_anno -> pirate ->  get_ags -> {gtdb uhgg}
panX -> subset_panx -> bakta_prot -> {gtdb uhgg}
}")



