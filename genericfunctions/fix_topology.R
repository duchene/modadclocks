
fix_tree_topology <- function(xml_path, input_tree, output_name){

  if(missing(input_tree)){
	input_tree <- '(Chelonia_mydas_5,Chelonia_mydas_3,((((Chelonia_mydas_1,Chelonia_mydas_6),(((((((Caretta_caretta_1,Caretta_caretta_4),Caretta_caretta_3),Caretta_caretta_2),((Lepidochelys_kempii_1,Lepidochelys_kempii_2),((Lepidochelys_olivacea_1,Lepidochelys_olivacea_3),Lepidochelys_olivacea_2))),((Eretmochelys_imbricata_2,Eretmochelys_imbricata_3),Eretmochelys_imbricata_1)),(Dermochelys_coriacea_4,(Dermochelys_coriacea_1,(Dermochelys_coriacea_2,Dermochelys_coriacea_3)))),(Natator_depressa_1,Natator_depressa_2))),(Chelonia_mydas_2,Chelonia_mydas_4)),Chelonia_mydas_7));'
  }

  if(missing(xml_path)){
	raw_lines <- readLines('turtles_strict.xml')
  }else{
	raw_lines <- readLines(xml_path)
	print('read xml file') 
  }

  f_name = gsub('[.]xml', '', xml_path)

  rand_start_lines <- grep('(<init id=\"RandomTree)|</init>', raw_lines)
  raw_lines[rand_start_lines[1]:rand_start_lines[2]] <- ''
  tree_start_lines <- grep('(<tree id=)|(</tree>)', raw_lines)
  raw_lines[tree_start_lines[1]:tree_start_lines[2]] <- ''
  raw_lines[tree_start_lines[1]] <- paste0('<input name=\"stateNode\" idref=\"Tree.t:', f_name, '\"/>')

  fix_tree_location <- grep('<run id=', raw_lines)

  tree_block <- paste0('<tree id=\"Tree.t:', f_name,'\" spec=\"beast.util.TreeParser\" newick=\"', input_tree, '\" taxa=\"@', f_name, '\" IsLabelledNewick=\"true\"/>')

  raw_lines <- c(raw_lines[1:(fix_tree_location-1)], tree_block, raw_lines[fix_tree_location:length(raw_lines)])

  # Remove operators
  ops_lines <- grep('WilsonBalding|SubtreeSlide|narrow[.]t|wide[.]t', raw_lines)
  raw_lines[ops_lines] <- ''


  if(missing(output_name)){
	output_name <- 'test_insert_tree.xml'
  }

  cat(raw_lines, sep = '\n', file = output_name)
  cat('wrote output to:', output_name, '\n')

}