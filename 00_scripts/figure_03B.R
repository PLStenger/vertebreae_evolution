# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005975","carbohydrate metabolic process",5.177832504046095,0.9572038772201659,0.09992703,"carbohydrate metabolic process"),
                     c("GO:0006629","lipid metabolic process",5.980051423570241,0.9565158708717862,0.12789186,"carbohydrate metabolic process"),
                     c("GO:0006099","tricarboxylic acid cycle",0.4936760212682554,0.9341599740463774,0.0832772,"tricarboxylic acid cycle"),
                     c("GO:0005978","glycogen biosynthetic process",0.08594552208766443,0.907289650794849,0.66692831,"tricarboxylic acid cycle"),
                     c("GO:0006094","gluconeogenesis",0.25553027855475857,0.9005981495250199,0.67250444,"tricarboxylic acid cycle"),
                     c("GO:0006509","membrane protein ectodomain proteolysis",0.020768896138831334,0.9477742153738016,0.03430977,"membrane protein ectodomain proteolysis"),
                     c("GO:0006260","DNA replication",1.4421509940592292,0.9117034334134572,0.1151447,"membrane protein ectodomain proteolysis"),
                     c("GO:0006285","base-excision repair, AP site formation",0.10580371763733157,0.8646418695982601,0.49165821,"membrane protein ectodomain proteolysis"),
                     c("GO:0006310","DNA recombination",1.827734239452524,0.9096835256835772,0.66164169,"membrane protein ectodomain proteolysis"),
                     c("GO:0006353","DNA-templated transcription termination",0.1805525451841636,0.8872053798511512,0.54875679,"membrane protein ectodomain proteolysis"),
                     c("GO:0006364","rRNA processing",1.5176702250774354,0.7863705663068109,0.69802204,"membrane protein ectodomain proteolysis"),
                     c("GO:0006367","transcription initiation at RNA polymerase II promoter",0.1920519861370695,0.886672296263961,0.60642105,"membrane protein ectodomain proteolysis"),
                     c("GO:0006396","RNA processing",4.402491558667012,0.85307499353138,0.47697115,"membrane protein ectodomain proteolysis"),
                     c("GO:0006397","mRNA processing",1.2042316958165515,0.8647092646294445,0.67661264,"membrane protein ectodomain proteolysis"),
                     c("GO:0006398","mRNA 3'-end processing by stem-loop binding and cleavage",0.019993569961569913,0.9006045901956001,0.2730445,"membrane protein ectodomain proteolysis"),
                     c("GO:0006457","protein folding",1.1925033490398829,0.8878624684467661,0.32158831,"membrane protein ectodomain proteolysis"),
                     c("GO:0006508","proteolysis",5.329242630377004,0.9260187888632558,0.39133953,"membrane protein ectodomain proteolysis"),
                     c("GO:0006511","ubiquitin-dependent protein catabolic process",1.2078769595261523,0.8933955459227878,0.43048359,"membrane protein ectodomain proteolysis"),
                     c("GO:0008380","RNA splicing",0.9583080778010037,0.8706919494898192,0.65672306,"membrane protein ectodomain proteolysis"),
                     c("GO:0009100","glycoprotein metabolic process",0.9596544378612005,0.9155447634790511,0.31348101,"membrane protein ectodomain proteolysis"),
                     c("GO:0009303","rRNA transcription",0.051247829640508086,0.8957521322859617,0.31709385,"membrane protein ectodomain proteolysis"),
                     c("GO:0016070","RNA metabolic process",8.105171248384922,0.8998690362778108,0.6180467,"membrane protein ectodomain proteolysis"),
                     c("GO:0022417","protein maturation by protein folding",0.022570606493610253,0.917039110525155,0.57970253,"membrane protein ectodomain proteolysis"),
                     c("GO:0035523","protein K29-linked deubiquitination",0.004469816945735684,0.9470689994331746,0.5070786,"membrane protein ectodomain proteolysis"),
                     c("GO:0043687","post-translational protein modification",1.8860683042179072,0.9261007675874033,0.21696733,"membrane protein ectodomain proteolysis"),
                     c("GO:0050650","chondroitin sulfate proteoglycan biosynthetic process",0.04048925592365198,0.8938213128228207,0.40690582,"membrane protein ectodomain proteolysis"),
                     c("GO:0051156","glucose 6-phosphate metabolic process",0.2974889621857917,0.9048547426791963,0.484459,"membrane protein ectodomain proteolysis"),
                     c("GO:0070646","protein modification by small protein removal",0.39107698516360445,0.9305108841925455,0.48050258,"membrane protein ectodomain proteolysis"),
                     c("GO:0071947","protein deubiquitination involved in ubiquitin-dependent protein catabolic process",0.014086322896720991,0.9114073053714167,0.6475366,"membrane protein ectodomain proteolysis"),
                     c("GO:1990168","protein K33-linked deubiquitination",0.004487046416341493,0.947057107804904,0.56976397,"membrane protein ectodomain proteolysis"),
                     c("GO:0006915","apoptotic process",0.3770497347375319,0.9740500444109896,0.01442917,"apoptotic process"),
                     c("GO:0008631","intrinsic apoptotic signaling pathway in response to oxidative stress",0.0038397105921517986,0.8761665363144386,0.69219823,"apoptotic process"),
                     c("GO:0043653","mitochondrial fragmentation involved in apoptotic process",0.003984930415829335,0.9061108055297008,0.6938338,"apoptotic process"),
                     c("GO:0007018","microtubule-based movement",0.5954751176662081,0.9599183965507159,0.01509482,"microtubule-based movement"),
                     c("GO:0051012","microtubule sliding",0.005003930534515774,0.969840628444691,0.56498176,"microtubule-based movement"),
                     c("GO:0051228","mitotic spindle disassembly",0.012639047365833005,0.8893306997134707,0.60210976,"microtubule-based movement"),
                     c("GO:0051256","mitotic spindle midzone assembly",0.015853574310288294,0.8809994152197264,0.69409251,"microtubule-based movement"),
                     c("GO:0007154","cell communication",9.236864411598072,0.9805157330685438,0.02242922,"cell communication"),
                     c("GO:0007155","cell adhesion",1.0284049482842672,0.9843198653023745,0.01597606,"cell adhesion"),
                     c("GO:0007528","neuromuscular junction development",0.02862307338213671,0.9250239748700347,0.01155453,"neuromuscular junction development"),
                     c("GO:0000002","mitochondrial genome maintenance",0.06554336753744247,0.9005698079271863,0.41231643,"neuromuscular junction development"),
                     c("GO:0000266","mitochondrial fission",0.03339317538700222,0.901199926224671,0.65271177,"neuromuscular junction development"),
                     c("GO:0006325","chromatin organization",1.6774440287109929,0.9031676765314288,0.36106458,"neuromuscular junction development"),
                     c("GO:0007005","mitochondrion organization",0.6594013763196482,0.8926404625634629,0.40277433,"neuromuscular junction development"),
                     c("GO:0007010","cytoskeleton organization",1.8978237258769566,0.8832923875595957,0.56675065,"neuromuscular junction development"),
                     c("GO:0016559","peroxisome fission",0.02904888744139457,0.9098321011701817,0.28289134,"neuromuscular junction development"),
                     c("GO:0030198","extracellular matrix organization",0.23021526352893723,0.9168997182810024,0.25779292,"neuromuscular junction development"),
                     c("GO:0032060","bleb assembly",0.0011543745305892267,0.9188692298429769,0.18788686,"neuromuscular junction development"),
                     c("GO:0033962","P-body assembly",0.043898229750658546,0.895755914453354,0.38353708,"neuromuscular junction development"),
                     c("GO:0045109","intermediate filament organization",0.06798502965758002,0.8916973624556741,0.41354099,"neuromuscular junction development"),
                     c("GO:0048268","clathrin coat assembly",0.03653386174314689,0.9097925356758892,0.31117528,"neuromuscular junction development"),
                     c("GO:0048312","intracellular distribution of mitochondria",0.006805640889294695,0.8791560733943872,0.58946724,"neuromuscular junction development"),
                     c("GO:0071788","endoplasmic reticulum tubular network maintenance",0.005227913652391295,0.9187056342614177,0.20363779,"neuromuscular junction development"),
                     c("GO:1905515","non-motile cilium assembly",0.05417437829055198,0.8837424928147185,0.57332702,"neuromuscular junction development"),
                     c("GO:0007586","digestion",0.03046170403107094,0.9136290656954805,-0,"digestion"),
                     c("GO:0001503","ossification",0.058122388412226014,0.9100905935479313,0.41288749,"digestion"),
                     c("GO:0001570","vasculogenesis",0.016510755546252736,0.8296391675809999,0.43429718,"digestion"),
                     c("GO:0006936","muscle contraction",0.11758129147287409,0.9014389196527368,0.65450471,"digestion"),
                     c("GO:0007601","visual perception",0.10847920828711939,0.9019568042535839,0.61012035,"digestion"),
                     c("GO:0016203","muscle attachment",0.003920935239293472,0.855357109616714,0.42527118,"digestion"),
                     c("GO:0021515","cell differentiation in spinal cord",0.012191081130081962,0.837735137354156,0.50190196,"digestion"),
                     c("GO:0021675","nerve development",0.021359620845316226,0.8486142602763164,0.50123956,"digestion"),
                     c("GO:0032835","glomerulus development",0.012124624600602412,0.8520283462505918,0.47133475,"digestion"),
                     c("GO:0043584","nose development",0.0015211161191985973,0.8669903807874583,0.33450301,"digestion"),
                     c("GO:0048484","enteric nervous system development",0.011912948247445324,0.8537109955190599,0.38284791,"digestion"),
                     c("GO:0060047","heart contraction",0.04380715969174213,0.9033438236575762,0.4227675,"digestion"),
                     c("GO:0060059","embryonic retina morphogenesis in camera-type eye",0.0012380605306745866,0.860095662878874,0.5614342,"digestion"),
                     c("GO:0071690","cardiac muscle myosin thick filament assembly",4.676570593005396E-05,0.7994502217654257,0.49990411,"digestion"),
                     c("GO:0009235","cobalamin metabolic process",0.24150056677574236,0.936616646649158,0.04209001,"cobalamin metabolic process"),
                     c("GO:0006011","UDP-glucose metabolic process",0.04875201775560941,0.87595511404449,0.42156521,"cobalamin metabolic process"),
                     c("GO:0006020","inositol metabolic process",0.10161203357423251,0.9356813258161137,0.25211295,"cobalamin metabolic process"),
                     c("GO:0006101","citrate metabolic process",0.0349167528591445,0.918004651187555,0.23285758,"cobalamin metabolic process"),
                     c("GO:0006108","malate metabolic process",0.0808972872001623,0.9117666010531721,0.30635422,"cobalamin metabolic process"),
                     c("GO:0006144","purine nucleobase metabolic process",0.2899621448839967,0.893772564163703,0.27438589,"cobalamin metabolic process"),
                     c("GO:0006171","cAMP biosynthetic process",0.05414730340817142,0.8569849404091864,0.57319204,"cobalamin metabolic process"),
                     c("GO:0006741","NADP biosynthetic process",0.04289399774963423,0.8680645178309327,0.48470802,"cobalamin metabolic process"),
                     c("GO:0019343","cysteine biosynthetic process via cystathionine",0.008747648361863777,0.8949797240803021,0.26697321,"cobalamin metabolic process"),
                     c("GO:0019674","NAD metabolic process",0.3681174849048916,0.8612781481466407,0.67381419,"cobalamin metabolic process"),
                     c("GO:0035338","long-chain fatty-acyl-CoA biosynthetic process",0.006621039418518167,0.8567819277253516,0.44022557,"cobalamin metabolic process"),
                     c("GO:0046037","GMP metabolic process",0.13243309513508175,0.8606502659850753,0.52852779,"cobalamin metabolic process"),
                     c("GO:0055129","L-proline biosynthetic process",0.09278069921228338,0.8849537412606929,0.45749557,"cobalamin metabolic process"),
                     c("GO:0015798","myo-inositol transport",0.022568145140666566,0.9499146674345645,-0,"myo-inositol transport"),
                     c("GO:0006626","protein targeting to mitochondrion",0.20164141720567424,0.8902098167033348,0.52110731,"myo-inositol transport"),
                     c("GO:0006811","monoatomic ion transport",4.917982551075166,0.9279546197709545,0.3621517,"myo-inositol transport"),
                     c("GO:0006836","neurotransmitter transport",0.14346487902868713,0.946471126655008,0.19075924,"myo-inositol transport"),
                     c("GO:0006865","amino acid transport",1.027831453048388,0.9374871970205372,0.24708658,"myo-inositol transport"),
                     c("GO:0006869","lipid transport",0.8830497501948285,0.9133035446055311,0.40882207,"myo-inositol transport"),
                     c("GO:0006891","intra-Golgi vesicle-mediated transport",0.13625065355074037,0.9078241737525351,0.20910405,"myo-inositol transport"),
                     c("GO:0006897","endocytosis",0.5813395677106135,0.8980572870955186,0.60588461,"myo-inositol transport"),
                     c("GO:0006900","vesicle budding from membrane",0.13789237596417964,0.8356582802626047,0.60650305,"myo-inositol transport"),
                     c("GO:0015677","copper ion import",0.038488175980434414,0.946953163475872,0.1691474,"myo-inositol transport"),
                     c("GO:0015730","propanoate transmembrane transport",0.001991234531442824,0.9277699954926564,0.59653222,"myo-inositol transport"),
                     c("GO:0016191","synaptic vesicle uncoating",0.007558814890062933,0.8473126238668576,0.63115813,"myo-inositol transport"),
                     c("GO:0016192","vesicle-mediated transport",2.591260690701909,0.9186271761083311,0.41239391,"myo-inositol transport"),
                     c("GO:0030705","cytoskeleton-dependent intracellular transport",0.26073850138360033,0.8991327236171948,0.68301253,"myo-inositol transport"),
                     c("GO:0032511","late endosome to vacuole transport via multivesicular body sorting pathway",0.11359389970410107,0.8879264702958427,0.62810709,"myo-inositol transport"),
                     c("GO:0034394","protein localization to cell surface",0.011829262247359966,0.9189750325751822,0.32452662,"myo-inositol transport"),
                     c("GO:0035494","SNARE complex disassembly",0.017921110782985417,0.8632935952102083,0.46742575,"myo-inositol transport"),
                     c("GO:0043252","sodium-independent organic anion transport",0.022760130670274156,0.9476993363873515,0.1637678,"myo-inositol transport"),
                     c("GO:0046907","intracellular transport",2.85589797514831,0.8792828369129714,0.55692439,"myo-inositol transport"),
                     c("GO:0071896","protein localization to adherens junction",0.0023013650023473923,0.9262739907738609,0.13718841,"myo-inositol transport"),
                     c("GO:0072583","clathrin-dependent endocytosis",0.07067774977797366,0.9107521692105054,0.69374477,"myo-inositol transport"),
                     c("GO:0072659","protein localization to plasma membrane",0.18653363283732313,0.897915016285612,0.51119726,"myo-inositol transport"),
                     c("GO:0097355","protein localization to heterochromatin",0.0040513869453088856,0.9197367778016211,0.29897635,"myo-inositol transport"),
                     c("GO:1990542","mitochondrial transmembrane transport",0.330613850101932,0.9098613266370995,0.26713411,"myo-inositol transport"),
                     c("GO:0016051","carbohydrate biosynthetic process",1.1871720585638565,0.9331657641108478,0.05722099,"carbohydrate biosynthetic process"),
                     c("GO:0006631","fatty acid metabolic process",1.7695330877461002,0.8465022205986088,0.62487897,"carbohydrate biosynthetic process"),
                     c("GO:0006694","steroid biosynthetic process",0.22341700669847356,0.9015479750850485,0.17029226,"carbohydrate biosynthetic process"),
                     c("GO:0016042","lipid catabolic process",1.1775407844952093,0.8806260127409732,0.50225527,"carbohydrate biosynthetic process"),
                     c("GO:0019367","fatty acid elongation, saturated fatty acid",0.03544348238909353,0.8622305363684543,0.50760125,"carbohydrate biosynthetic process"),
                     c("GO:0019369","arachidonate metabolic process",0.04666725181230648,0.8764540768103014,0.38454233,"carbohydrate biosynthetic process"),
                     c("GO:0019432","triglyceride biosynthetic process",0.0858569133816917,0.889929492176428,0.49222147,"carbohydrate biosynthetic process"),
                     c("GO:0030148","sphingolipid biosynthetic process",0.24209375283517093,0.8796344137493098,0.53429597,"carbohydrate biosynthetic process"),
                     c("GO:0030497","fatty acid elongation",0.0736535254868913,0.8586674001389712,0.6473653,"carbohydrate biosynthetic process"),
                     c("GO:0042761","very long-chain fatty acid biosynthetic process",0.058110081647507576,0.860916340996123,0.61660944,"carbohydrate biosynthetic process"),
                     c("GO:0016477","cell migration",0.48278207313949656,0.9785980337700029,0.01478176,"cell migration"),
                     c("GO:1904969","slow muscle cell migration",3.692029415530576E-05,0.9865638956263931,0.47080365,"cell migration"),
                     c("GO:0019254","carnitine metabolic process, CoA-linked",0.0035615777095151626,0.9604899931847526,0.03028627,"carnitine metabolic process, CoA-linked"),
                     c("GO:0046314","phosphocreatine biosynthetic process",0.014061709367284122,0.9205105259659444,0.46805621,"carnitine metabolic process, CoA-linked"),
                     c("GO:0030048","actin filament-based movement",0.08138217373006863,0.9623541831342987,0.01256947,"actin filament-based movement"),
                     c("GO:0043149","stress fiber assembly",0.006406901712417393,0.8928244831284233,0.57058383,"actin filament-based movement"),
                     c("GO:0032259","methylation",2.6585910004764686,0.9737125464371951,0.0540798,"methylation"),
                     c("GO:0033627","cell adhesion mediated by integrin",0.05005899616870724,0.9709802583435081,0.01207616,"cell adhesion mediated by integrin"),
                     c("GO:0007157","heterophilic cell-cell adhesion via plasma membrane cell adhesion molecules",0.036339414860595616,0.9700866952610023,0.58378397,"cell adhesion mediated by integrin"),
                     c("GO:0007160","cell-matrix adhesion",0.09374554956620869,0.9698152080536864,0.62174985,"cell adhesion mediated by integrin"),
                     c("GO:0042554","superoxide anion generation",0.0302795639132381,0.957309707357306,0.06089005,"superoxide anion generation"),
                     c("GO:0042593","glucose homeostasis",0.07836947772699569,0.9941614978793794,-0,"glucose homeostasis"),
                     c("GO:0007035","vacuolar acidification",0.05759565888227698,0.9064431071671536,0.56691239,"glucose homeostasis"),
                     c("GO:0045197","establishment or maintenance of epithelial cell apical/basal polarity",0.051956699288289956,0.9833413616099265,0.01211256,"establishment or maintenance of epithelial cell apical/basal polarity"),
                     c("GO:0046185","aldehyde catabolic process",0.1262994035994136,0.913280623281827,0.09320059,"aldehyde catabolic process"),
                     c("GO:0006538","glutamate catabolic process",0.04895138734404807,0.876328219367958,0.63763456,"aldehyde catabolic process"),
                     c("GO:0006572","tyrosine catabolic process",0.04589930969387612,0.8776797976160199,0.60750873,"aldehyde catabolic process"),
                     c("GO:0006574","valine catabolic process",0.060145620531936776,0.8778286456195115,0.46638159,"aldehyde catabolic process"),
                     c("GO:0009395","phospholipid catabolic process",0.11020707805358769,0.8507768905185186,0.48732639,"aldehyde catabolic process"),
                     c("GO:0034638","phosphatidylcholine catabolic process",0.0023604374729958816,0.8808145502647855,0.54215157,"aldehyde catabolic process"),
                     c("GO:0061709","reticulophagy",0.04193160874865259,0.9129424149525055,0.45473981,"aldehyde catabolic process"),
                     c("GO:0048066","developmental pigmentation",0.023939118730300253,1,-0,"developmental pigmentation"),
                     c("GO:0048511","rhythmic process",0.13974577473077598,1,0,"rhythmic process"),
                     c("GO:0048870","cell motility",0.8376747086779577,0.9845996546371083,-0,"cell motility"),
                     c("GO:0050795","regulation of behavior",0.015243158780253906,0.926148472452721,0.08777637,"regulation of behavior"),
                     c("GO:0000381","regulation of alternative mRNA splicing, via spliceosome",0.07535185901803537,0.9114225729471794,0.21165947,"regulation of behavior"),
                     c("GO:0006357","regulation of transcription by RNA polymerase II",4.1876991326832185,0.8761129092626275,0.36822218,"regulation of behavior"),
                     c("GO:0007166","cell surface receptor signaling pathway",1.7427092633657988,0.827793999900378,0.53421354,"regulation of behavior"),
                     c("GO:0007179","transforming growth factor beta receptor signaling pathway",0.05119614122869066,0.8605040479417645,0.5121069,"regulation of behavior"),
                     c("GO:0007186","G protein-coupled receptor signaling pathway",1.464221945905271,0.8304196916160612,0.36267237,"regulation of behavior"),
                     c("GO:0007188","adenylate cyclase-modulating G protein-coupled receptor signaling pathway",0.34707291723636735,0.8447754174480238,0.44702642,"regulation of behavior"),
                     c("GO:0007214","gamma-aminobutyric acid signaling pathway",0.030321406913280776,0.8896501862828917,0.1422804,"regulation of behavior"),
                     c("GO:0007219","Notch signaling pathway",0.0884413339725631,0.8552556346183009,0.56469551,"regulation of behavior"),
                     c("GO:0007229","integrin-mediated signaling pathway",0.08566246649914042,0.8555730309088536,0.54360877,"regulation of behavior"),
                     c("GO:0007249","canonical NF-kappaB signal transduction",0.028052039499201313,0.8671815343882944,0.55664143,"regulation of behavior"),
                     c("GO:0008277","regulation of G protein-coupled receptor signaling pathway",0.06807856106944013,0.8774349001578886,0.52211654,"regulation of behavior"),
                     c("GO:0008285","negative regulation of cell population proliferation",0.12504657495107693,0.8893741890974719,0.17177237,"regulation of behavior"),
                     c("GO:0010629","negative regulation of gene expression",1.0018051070218412,0.8608083774029207,0.53324106,"regulation of behavior"),
                     c("GO:0010634","positive regulation of epithelial cell migration",0.010221998775132322,0.8908895712152891,0.66701115,"regulation of behavior"),
                     c("GO:0010739","positive regulation of protein kinase A signaling",0.002293980943516331,0.8771442069382195,0.44823279,"regulation of behavior"),
                     c("GO:0010747","positive regulation of long-chain fatty acid import across plasma membrane",0.0018164784724410432,0.8989091596000276,0.56496602,"regulation of behavior"),
                     c("GO:0010811","positive regulation of cell-substrate adhesion",0.026836131145019912,0.8985874742514867,0.34402406,"regulation of behavior"),
                     c("GO:0014911","positive regulation of smooth muscle cell migration",0.0034508168270492453,0.9034770051625172,0.10200035,"regulation of behavior"),
                     c("GO:0016055","Wnt signaling pathway",0.19157448366599422,0.847119008333302,0.59889884,"regulation of behavior"),
                     c("GO:0019722","calcium-mediated signaling",0.1524365105084264,0.8511427776473897,0.65500853,"regulation of behavior"),
                     c("GO:0030177","positive regulation of Wnt signaling pathway",0.06665343771504534,0.851805807967032,0.63494124,"regulation of behavior"),
                     c("GO:0030336","negative regulation of cell migration",0.04620944016478069,0.8873325581639285,0.69764339,"regulation of behavior"),
                     c("GO:0031087","deadenylation-independent decapping of nuclear-transcribed mRNA",0.020741821256450774,0.802324822366994,0.55809235,"regulation of behavior"),
                     c("GO:0031507","heterochromatin formation",0.11730561994318113,0.8110196216249623,0.64327644,"regulation of behavior"),
                     c("GO:0032436","positive regulation of proteasomal ubiquitin-dependent protein catabolic process",0.11581896276519417,0.8870553898307934,0.41670478,"regulation of behavior"),
                     c("GO:0032483","regulation of Rab protein signal transduction",0.0195012993728325,0.8811693106220312,0.50057912,"regulation of behavior"),
                     c("GO:0033173","calcineurin-NFAT signaling cascade",0.011418216305764227,0.8743653798226252,0.25244228,"regulation of behavior"),
                     c("GO:0034587","piRNA processing",0.01849952872475187,0.8061234125018716,0.55325473,"regulation of behavior"),
                     c("GO:0035024","negative regulation of Rho protein signal transduction",0.00890517495025975,0.8645998836832135,0.66897802,"regulation of behavior"),
                     c("GO:0035329","hippo signaling",0.028901206264773348,0.8706206794218854,0.34710849,"regulation of behavior"),
                     c("GO:0035556","intracellular signal transduction",4.385343312708344,0.8123984423776763,0.61586777,"regulation of behavior"),
                     c("GO:0035567","non-canonical Wnt signaling pathway",0.03480845332962227,0.8624260637686214,0.26995862,"regulation of behavior"),
                     c("GO:0038093","Fc receptor signaling pathway",0.004715952240104389,0.8492495090573131,0.44397717,"regulation of behavior"),
                     c("GO:0040036","regulation of fibroblast growth factor receptor signaling pathway",0.017460837782515937,0.8840636051846666,0.4339994,"regulation of behavior"),
                     c("GO:0042058","regulation of epidermal growth factor receptor signaling pathway",0.03246278397428851,0.883090289493336,0.47953872,"regulation of behavior"),
                     c("GO:0042127","regulation of cell population proliferation",0.401650957409684,0.9100507882757504,0.1359632,"regulation of behavior"),
                     c("GO:0042981","regulation of apoptotic process",0.49277516609086597,0.8986877268159569,0.20266419,"regulation of behavior"),
                     c("GO:0043408","regulation of MAPK cascade",0.2982839791866026,0.86096358740337,0.63895379,"regulation of behavior"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.531020603443994,0.8526544591637082,0.64398089,"regulation of behavior"),
                     c("GO:0046427","positive regulation of receptor signaling pathway via JAK-STAT",0.014332458191089695,0.8665619722240576,0.52542138,"regulation of behavior"),
                     c("GO:0048015","phosphatidylinositol-mediated signaling",0.08115080655336206,0.86140810368451,0.3916997,"regulation of behavior"),
                     c("GO:0048167","regulation of synaptic plasticity",0.061676582062910115,0.8670579050329924,0.52549383,"regulation of behavior"),
                     c("GO:0050709","negative regulation of protein secretion",0.009719882774620164,0.8902906353015777,0.52378643,"regulation of behavior"),
                     c("GO:0050793","regulation of developmental process",1.4293716495756061,0.9062992582253933,0.19842072,"regulation of behavior"),
                     c("GO:0050810","regulation of steroid biosynthetic process",0.016104632310544372,0.9250301833599374,0.11057902,"regulation of behavior"),
                     c("GO:0051209","release of sequestered calcium ion into cytosol",0.07223086348544018,0.819230956222316,0.56599593,"regulation of behavior"),
                     c("GO:0051726","regulation of cell cycle",0.9785034787039558,0.9032450524538723,0.21832891,"regulation of behavior"),
                     c("GO:0051781","positive regulation of cell division",0.05264095540663495,0.894143759454595,0.40345783,"regulation of behavior"),
                     c("GO:0061035","regulation of cartilage development",0.014509675603035163,0.9210054080879893,0.53498614,"regulation of behavior"),
                     c("GO:0090036","regulation of protein kinase C signaling",0.005636498241043346,0.8918233336543853,0.13697758,"regulation of behavior"),
                     c("GO:1902514","regulation of calcium ion transmembrane transport via high voltage-gated calcium channel",0.003453278179992932,0.9182020688039079,0.10164689,"regulation of behavior"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.0070099331836207205,0.9023907884979285,0.49646475,"regulation of behavior"),
                     c("GO:0051301","cell division",1.377833380287743,0.9839024644434342,0.0168732,"cell division"),
                     c("GO:0055018","regulation of cardiac muscle fiber development",0.0001698333531144065,0.9388063464672516,-0,"regulation of cardiac muscle fiber development"),
                     c("GO:0045661","regulation of myoblast differentiation",0.015134859250731674,0.9225232576545309,0.50055338,"regulation of cardiac muscle fiber development"),
                     c("GO:0070085","glycosylation",0.8286218525510767,0.9765120826540802,-0,"glycosylation"),
                     c("GO:0070121","Kupffer's vesicle development",0.001617108884002392,0.898249982701115,-0,"Kupffer's vesicle development"),
                     c("GO:0000902","cell morphogenesis",0.5643833072815533,0.8409761557733669,0.68788457,"Kupffer's vesicle development"),
                     c("GO:0002009","morphogenesis of an epithelium",0.09926636421889876,0.850341689314833,0.55247651,"Kupffer's vesicle development"),
                     c("GO:0007548","sex differentiation",0.055557658644904105,0.8751568511791815,0.33864486,"Kupffer's vesicle development"),
                     c("GO:0009888","tissue development",0.631041667702486,0.8520684559426898,0.49950736,"Kupffer's vesicle development"),
                     c("GO:0030316","osteoclast differentiation",0.01474596548562912,0.8538175828486763,0.384517,"Kupffer's vesicle development"),
                     c("GO:0030719","P granule organization",0.0036157274742762775,0.7892686561218759,0.61892027,"Kupffer's vesicle development"),
                     c("GO:0048513","animal organ development",1.042075302533505,0.8462007933147895,0.68737042,"Kupffer's vesicle development"),
                     c("GO:0048813","dendrite morphogenesis",0.01966867137300322,0.7774884483238135,0.60167101,"Kupffer's vesicle development"),
                     c("GO:0051216","cartilage development",0.04374562586814995,0.8362830008856957,0.68719552,"Kupffer's vesicle development"),
                     c("GO:0070593","dendrite self-avoidance",0.01952837425521306,0.8244259202798864,0.51899261,"Kupffer's vesicle development"),
                     c("GO:0070602","regulation of centromeric sister chromatid cohesion",0.002087227296246619,0.920791659574192,0.09903122,"regulation of centromeric sister chromatid cohesion"),
                     c("GO:0007026","negative regulation of microtubule depolymerization",0.015287463133240272,0.8868277432856557,0.5305296,"regulation of centromeric sister chromatid cohesion"),
                     c("GO:0045724","positive regulation of cilium assembly",0.009185769185840072,0.8888070348260727,0.64965176,"regulation of centromeric sister chromatid cohesion"),
                     c("GO:0090141","positive regulation of mitochondrial fission",0.01299348218972394,0.881975196821913,0.68028989,"regulation of centromeric sister chromatid cohesion"),
                     c("GO:1903358","regulation of Golgi organization",0.012395373424407987,0.9125538055030729,0.47710294,"regulation of centromeric sister chromatid cohesion"),
                     c("GO:1903373","positive regulation of endoplasmic reticulum tubular network organization",0.005496201123253184,0.8924565792315443,0.50166843,"regulation of centromeric sister chromatid cohesion"),
                     c("GO:0071356","cellular response to tumor necrosis factor",0.04704630016563428,0.9435719773184009,-0,"cellular response to tumor necrosis factor"),
                     c("GO:0006805","xenobiotic metabolic process",0.14736120073854372,0.8875956229807578,0.61252549,"cellular response to tumor necrosis factor"),
                     c("GO:0006935","chemotaxis",0.48949664396987486,0.9308986674743247,0.46382379,"cellular response to tumor necrosis factor"),
                     c("GO:0006952","defense response",1.1245700077941203,0.9285286458787003,0.54266122,"cellular response to tumor necrosis factor"),
                     c("GO:0006955","immune response",0.6332667307635791,0.9334327755584282,0.27394137,"cellular response to tumor necrosis factor"),
                     c("GO:0006979","response to oxidative stress",0.7667902052527142,0.9308746374484016,0.48029419,"cellular response to tumor necrosis factor"),
                     c("GO:0007212","G protein-coupled dopamine receptor signaling pathway",0.03525149685948594,0.8523901478959445,0.68094313,"cellular response to tumor necrosis factor"),
                     c("GO:0009408","response to heat",0.3335994712226244,0.935492485642144,0.22065342,"cellular response to tumor necrosis factor"),
                     c("GO:0034599","cellular response to oxidative stress",0.33220388410355384,0.9055982794539371,0.57145744,"cellular response to tumor necrosis factor"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.2589072547934972,0.9178399630979546,0.43556853,"cellular response to tumor necrosis factor"),
                     c("GO:0071277","cellular response to calcium ion",0.06120154094477852,0.9250116803727395,0.4734166,"cellular response to tumor necrosis factor"),
                     c("GO:1904613","cellular response to 2,3,7,8-tetrachlorodibenzodioxine",0.0014005098249579317,0.9389464530917708,0.31962482,"cellular response to tumor necrosis factor"),
                     c("GO:1900825","regulation of membrane depolarization during cardiac muscle cell action potential",0.001892780413695342,0.9215119903403777,0.09853873,"regulation of membrane depolarization during cardiac muscle cell action potential"),
                     c("GO:0007204","positive regulation of cytosolic calcium ion concentration",0.08884745720827146,0.9093945917020109,0.3788827,"regulation of membrane depolarization during cardiac muscle cell action potential"),
                     c("GO:0007596","blood coagulation",0.07478574784098735,0.7934043829352645,0.47681318,"regulation of membrane depolarization during cardiac muscle cell action potential"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "uniqueness",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()


# Custom from now:

###############################################################################
# 1. Packages and settings
###############################################################################

library(tidyverse)
library(ggrepel)

# install.packages(c("tidyverse", "ggrepel"))

input_script <- "Revigo_BP_TreeMap.R"
output_dir <- "REVIGO_alternative_figures"

dir.create(output_dir, showWarnings = FALSE)

n_groups <- 15       # Groups displayed in overview figures
n_detail_groups <- 6 # Groups displayed in the detailed figure
n_terms <- 6         # Terms displayed per detailed group

###############################################################################
# 2. Extract data without running the original plotting commands
###############################################################################

script_expressions <- parse(input_script)
data_environment <- new.env(parent = baseenv())

for (expression in script_expressions) {
  if (
    is.call(expression) &&
    as.character(expression[[1]]) %in% c("<-", "=") &&
    is.symbol(expression[[2]]) &&
    as.character(expression[[2]]) %in% c("revigo.names", "revigo.data")
  ) {
    eval(expression, envir = data_environment)
  }
}

go <- as_tibble(
  as.data.frame(
    data_environment$revigo.data,
    stringsAsFactors = FALSE
  ),
  .name_repair = "minimal"
)

names(go) <- data_environment$revigo.names


go <- go %>%
  mutate(
    across(c(frequency, uniqueness, dispensability), as.numeric)
  ) %>%
  distinct(term_ID, .keep_all = TRUE)

stopifnot(
  all(is.finite(go$frequency)),
  all(is.finite(go$uniqueness)),
  all(is.finite(go$dispensability))
)

###############################################################################
# 3. Summarize REVIGO groups
###############################################################################

groups <- go %>%
  count(representative, name = "n_terms") %>%
  arrange(desc(n_terms), representative)

top_groups <- groups %>%
  slice_head(n = n_groups)

# Use identical group ordering across overview figures.
group_order <- rev(top_groups$representative)

# Retrieve the actual representative term, not an average of its group.
representatives <- go %>%
  filter(description == representative) %>%
  inner_join(top_groups, by = "representative") %>%
  mutate(
    representative = factor(representative, levels = group_order)
  )

write_tsv(go, file.path(output_dir, "REVIGO_terms.tsv"))
write_tsv(groups, file.path(output_dir, "REVIGO_group_summary.tsv"))

###############################################################################
# 4. Shared publication-style theme and export function
###############################################################################

theme_publication <- function() {
  theme_classic(base_size = 10, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      axis.text = element_text(colour = "grey20"),
      axis.title = element_text(size = 10),
      axis.ticks = element_line(linewidth = 0.3),
      axis.line = element_line(linewidth = 0.3),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0),
      plot.margin = margin(10, 18, 10, 10)
    )
}

save_figure <- function(plot, name, width, height) {

  ggsave(
    file.path(output_dir, paste0(name, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = grDevices::pdf,
    useDingbats = FALSE,
    bg = "white"
  )

  ggsave(
    file.path(output_dir, paste0(name, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
}

###############################################################################
# 5. Figure 1: group-size lollipop plot
###############################################################################

p1 <- top_groups %>%
  mutate(
    representative = factor(representative, levels = group_order)
  ) %>%
  ggplot(aes(x = n_terms, y = representative)) +
  geom_segment(
    aes(x = 0, xend = n_terms, yend = representative),
    linewidth = 0.7,
    colour = "#BDD5DE"
  ) +
  geom_point(size = 3.2, colour = "#236B80") +
  geom_text(
    aes(label = n_terms),
    nudge_x = max(top_groups$n_terms) * 0.035,
    hjust = 0,
    size = 3
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, 43)) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.13)),
    breaks = scales::breaks_pretty(n = 5)
  ) +
  labs(
    title = "Functional composition of REVIGO groups",
    subtitle = paste("The", nrow(top_groups), "largest groups"),
    x = "Number of GO terms in the supplied REVIGO output",
    y = NULL
  ) +
  theme_publication()

save_figure(p1, "01_group_lollipop", width = 9, height = 7)

###############################################################################
# 6. Figure 2: representative-term dot plot
###############################################################################

p2 <- ggplot(
  representatives,
  aes(x = uniqueness, y = representative)
) +
  geom_point(
    aes(size = n_terms, colour = frequency),
    alpha = 0.95
  ) +
  scale_size_area(
    max_size = 9,
    name = "GO terms\nin group"
  ) +
  scale_colour_viridis_c(
    option = "C",
    end = 0.9,
    name = "REVIGO\nfrequency (%)"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, 43)) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  labs(
    title = "Semantic uniqueness of representative GO terms",
    subtitle = "Groups ordered by the number of associated GO terms",
    x = "Representative-term uniqueness",
    y = NULL
  ) +
  theme_publication() +
  theme(panel.grid.major.x = element_line(colour = "grey92"))

save_figure(p2, "02_representative_dotplot", width = 10, height = 7)

###############################################################################
# 7. Figure 3: detailed dot plots within the largest groups
###############################################################################

detail_groups <- groups %>%
  slice_head(n = n_detail_groups) %>%
  pull(representative)

detail <- go %>%
  filter(representative %in% detail_groups) %>%
  group_by(representative) %>%
  # Keep the least dispensable terms within each selected group.
  slice_min(
    order_by = dispensability,
    n = n_terms,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    representative = factor(
      representative,
      levels = detail_groups
    )
  ) %>%
  arrange(representative, desc(dispensability)) %>%
  mutate(
    term_label = paste0(description, " [", term_ID, "]"),
    term_label = factor(term_label, levels = unique(term_label))
  )

p3 <- ggplot(
  detail,
  aes(x = dispensability, y = term_label)
) +
  geom_segment(
    aes(x = 0, xend = dispensability, yend = term_label),
    colour = "grey85",
    linewidth = 0.5
  ) +
  geom_point(
    aes(size = frequency, colour = uniqueness),
    alpha = 0.95
  ) +
  facet_wrap(
    vars(representative),
    ncol = 2,
    scales = "free_y",
    labeller = label_wrap_gen(width = 35)
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, 43)) +
  scale_x_continuous(
    limits = c(0, max(0.7, detail$dispensability)),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  scale_size_area(
    max_size = 5,
    name = "REVIGO\nfrequency (%)"
  ) +
  scale_colour_viridis_c(
    option = "D",
    end = 0.85,
    name = "Uniqueness"
  ) +
  labs(
    title = "GO terms within the largest REVIGO groups",
    subtitle = "Up to six terms per group, selected by lowest dispensability",
    x = "Dispensability",
    y = NULL
  ) +
  theme_publication() +
  theme(
    axis.text.y = element_text(size = 7),
    panel.spacing = unit(1.4, "lines"),
    legend.position = "bottom"
  )

save_figure(p3, "03_group_details", width = 15, height = 10)

###############################################################################
# 8. Figure 4: uniqueness versus dispensability
###############################################################################

labelled_terms <- go %>%
  filter(description == representative) %>%
  inner_join(
    groups %>% slice_head(n = 10),
    by = "representative"
  )

p4 <- ggplot(go, aes(x = uniqueness, y = dispensability)) +
  geom_point(
    colour = "grey65",
    size = 1.8,
    alpha = 0.5
  ) +
  geom_point(
    data = labelled_terms,
    aes(size = n_terms),
    colour = "#236B80",
    alpha = 0.95
  ) +
  geom_text_repel(
    data = labelled_terms,
    aes(label = str_wrap(description, 28)),
    size = 2.8,
    seed = 123,
    box.padding = 0.6,
    point.padding = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    colour = "#17495A",
    segment.colour = "grey65"
  ) +
  scale_size_area(
    max_size = 7,
    name = "GO terms\nin group"
  ) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = c(0.20, 0.12))) +
  labs(
    title = "Gene Ontology terms of the duplicated gene copies following Ts3R",
    subtitle = "Semantic structure of Gene Ontology terms",
    x = "Uniqueness",
    y = "Dispensability"
  ) +
  theme_publication()

save_figure(p4, "04_semantic_overview", width = 10, height = 7)

###############################################################################
# 9. Display figures
###############################################################################

print(p1)
print(p2)
print(p3)
print(p4)

