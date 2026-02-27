(* ::Package:: *)

(* ::Title:: *)
(*Cell Related Functions*)


(* ::Input::Initialization:: *)
BeginPackage["cellRelatedFunctions`"]


(* ::Input::Initialization:: *)
Unprotect["`*"];ClearAll["`*"];
(*unprotect and clear out old definitionsn-allows repeated loading and alterations*)


(* ::Input::Initialization:: *)
cellRelatedFunctions::usage="This package is where all functions that were to capture a cell behavior, rule, or biological first principle such as cell cycle, quiencence, division and so on
Last Edit was May 2024 : 
added CNV+ploidy
and changed line 345 newCell[[7,2]]=Total[cellRelatedFunctions`copiesT[# (*activegenes*),rnapol,transcriptionTime,genomeLength]&/@newCell[[6,2]]];\[IndentingNewLine]to\[IndentingNewLine]newCell[[7,2]]=cellRelatedFunctions`copiesT[Total[newCell[[6,2]]](*activegenes*),rnapol,transcriptionTime,genomeLength];
";


(* ::Input::Initialization:: *)
cellcycleChangeOPT::usage="Function that changes a single phase in the cell cycle duration of a cell"


(* ::Input::Initialization:: *)
cellcycle01::usage="Random Real number between [0 and value] + value2 + value3"


(* ::Input::Initialization:: *)
cellcycle02::usage="Random Variate of Poison Distribution[ mean] + increment"


(* ::Input::Initialization:: *)
cellcycle03::usage="Random Variate of a Normal Distribution[ mean, sigma] + increment"


(* ::Input::Initialization:: *)
cellcycle04::usage="Random Variate of a Normal Distribution[ mean+increment*sigma, sigma] + increment"


(* ::Input::Initialization:: *)
transcriptionT::usage="calculates the total transcription time allocated in the cell cycle"


(* ::Input::Initialization:: *)
parentalTranscripts::usage="calculates the probaility and keeps a certain number of transcripts from parental cell"


(* ::Input::Initialization:: *)
parentalTranscripts1::usage="calculates the probaility and keeps a certain number of transcripts from parental cell"


(* ::Input::Initialization:: *)
parentalTranscripts2::usage="calculates the probaility and keeps a certain number of transcripts from parental cell"


(* ::Input::Initialization:: *)
parentalTranscriptsAssociation::usage="Takes a list of the transcripts (assumes the list holds every gene and the counts of the transcripts {3,1,1} represent 3 copies of gene1 represented as 1->3 2->1 and 3->1) then calculates the probaility of keeping a certain number of transcripts from the list and returns what remains
eg. parentalTranscriptsAssociation[{3,1,1},0.1] return {2,0,0}
parentalTranscriptsAssociation[{3,1,1},0.7] return {2,1,1}"


(* ::Input::Initialization:: *)
activeGenome::usage="a basic function that captures random effects of the grN and determines the active genome by using {{position,replacement rule, gene length}}"


(* ::Input::Initialization:: *)
competitionTranscripts::usage="
There are two versions of this
:=[ploidy_,grN_,genome_,transcriptList_]->assumes transcript list is a collection of transcripts {1,2,1,1,3,2,2,3}
transcripts compete to turn on or off a gene based on a grn, the function returns a list {{position,replacement rule, gene length}}
while 
:=[ploidycnv_,grN_,transcriptList_] assumes totals {3,3,2}
transcripts compete to turn on or off a gene based on a grn, the function returns a list of replacement rules given a grn eg {{1\[Rule]1,2\[Rule]1,3\[Rule]1},{},{2\[Rule]0,3\[Rule]1}}\[IndentingNewLine]we select tasks given the number of transcripts for each gene {3,1,1} results in {{2\[Rule]1,1\[Rule]1,1\[Rule]1},{},{2\[Rule]0}}\[IndentingNewLine]then we group them by gene\[IndentingNewLine]and select based on gene ploidy number\[IndentingNewLine]for 1 = {2\[Rule]1,1\[Rule]1}, for 2 ={2\[Rule]1,2\[Rule]0,1\[Rule]1,1\[Rule]1}
"


(* ::Input::Initialization:: *)
copiesT::usage="
calculates the expected number of transcript copies,  returns a list with the number of trancripts, there are two versions of this function
:=[ploidy_,cellactivegeneList_,rnaPolRate_,transcriptionTime_,genomeLength_]
 uses the ploidy parameter and assumes the genes are copied propotional to the ploidy
example
genome={{1,1},{2,2},{3,3}}
genomeLength={1,2,3}
ploidy=1
copiesT[ploidy,{1,1,1},{0,0},3,genome]
={3,1,1}
ploidy=2
copiesT[ploidy,{1,1,1},{0,0},3,genome]
={6,2,2}

:=[cellactivegeneList_,rnaPolRate_,transcriptionTime_,genomeLength_]
No ploidy parameter,will need to run it for each copy of the genome but accomodates different active gene copies calculates the expected number of transcript copies,  returns a list with the number of trancripts
example
genome={{1,1},{2,2},{3,3}}
genomeLength={1,2,3}
copiesT[{1,1,1},{0,0},3,genomeLength]
={3,1,1}
copiesT[{1,1,1},{0,0},1,genomeLength]
={1,0,0}
if rnapolRate[[2]]>0 mean that there is a reintiation for a given distance before the end of the gene is reached
copiesT[{1,1,1},{0,0.5},3,genomeLength]
={8,3,1}
"


(* ::Input::Initialization:: *)
transcriptsT::usage="calculates the total number of transcripts in a cell based on expected copies and active genome"


(* ::Input::Initialization:: *)
randDistofT::usage="divided the transcripts randomly between the progeny cells"


(* ::Input::Initialization:: *)
randDistofT2::usage="given a list of transcripts and their amount 
eg \[LeftAssociation]1\[Rule]3,2\[Rule]1,3\[Rule]1\[RightAssociation]
takes values eg {3,1,1} and 
returns a random division of the transcripts eg {{2,0,1},{1,1,0}}
"


(* ::Input::Initialization:: *)
getCellPhase::usage="for each cell phase is calculated that fits within the time"


(* ::Input::Initialization:: *)
oneCellTranscriptCompetition::usage="each transcript will compete to interact with the genome"


(* ::Input::Initialization:: *)
transcriptDegredation::usage="a proportion of the transcripts will be lost"


(* ::Input::Initialization:: *)
transcriptCopyNumber::usage="calculates the total number of transcripts in a cell based on expected copies and active genome within the cell cycle time"


(* ::Input::Initialization:: *)
Begin["`Private`"]


(* ::Input::Initialization:: *)
(*Context[zPt];
Context[qPt];*)


(* ::Input::Initialization:: *)
cellcycle01[val_]:=RandomReal[val[[1]],WorkingPrecision->2]+val[[2]]+val[[3]]


(* ::Input::Initialization:: *)
cellcycle02[val_]:=RandomVariate[PoissonDistribution[val[[1]]]]+val[[2]]+val[[3]]


(* ::Input::Initialization:: *)
cellcycle03[val_]:=RandomVariate[NormalDistribution[val[[1]],val[[2]]]]+val[[3]]


(* ::Input::Initialization:: *)
cellcycle04[val_]:=RandomVariate[NormalDistribution[val[[1]]+(val[[2]]*val[[3]]),val[[2]]]]+val[[3]]


(* ::Input::Initialization:: *)
cellcycleChangeOPT[val_,opt_]:=Which[opt==0,val[[3]]+val[[1]],
opt==-1,cellcycle01[val],
opt==-2,cellcycle02[val],
opt==-3, If[(val[[2]]>0),cellcycle03[val] ,val[[1]]](*Normal Distribution[ mean+(increment *sigma), sigma] + increment*),opt==-4,If[(val[[2]]>0),cellcycle04[val] ,val[[1]]],
opt==-5,(Max[{val[[1]],val[[2]]}]*val[[3]])+Max[{val[[1]],val[[2]]}],
opt==-6,If[(val[[2]]>0),{val[[1]]-#,#}&@cellcycle03[val],{(Abs[val[[3]]]*val[[1]]),(1-Abs[val[[3]]])*val[[1]]}],
opt==-7,(({1-#,#}&@cellcycle03[{val[[3]],val[[2]],0}]))*val[[1]],
opt==-8,(({1-#,cellcycle03[{#,val[[2]],0}]}&@cellcycle03[{val[[3]],val[[2]],0}]))*val[[1]]
(*{split between two daughter cells 0\[Rule]val,1-*),
opt=="stemlike",(If[val[[5]]==1,val[[4]],{val[[1]],val[[1]]}])
(*{the example tries to capture stem cells splitting assymetrically basically val[[4]] is a parameter you set for {stem cell cell cycle, non stem cell cell cycle}}*)
]


(* ::Input::Initialization:: *)
transcriptionT[cellcycledurations_,cyclePhase_]:=Accumulate[cellcycledurations/. _?Negative->0][[cyclePhase]]


(* ::Input::Initialization:: *)
parentalTranscripts1[transcripts_,probability_]:=RandomSample[transcripts,UpTo[((Floor[Length[transcripts]*RandomReal[{0,probability}]])/. _?Negative->0)]
] 


(* ::Input::Initialization:: *)
parentalTranscripts[transcripts_,probability_]:=RandomSample[transcripts,UpTo[((Floor[Length[transcripts]*probability])/. _?Negative->0)]
] 


(* ::Input::Initialization:: *)
parentalTranscripts2[transcripts_,probability_]:=RandomSample[transcripts,UpTo[((Floor[Length[transcripts]*If[RandomReal[]<=probability,probability,1]])/. _?Negative->0)]
] 


(* ::Input::Initialization:: *)
parentalTranscriptsAssociation[transcripts_,probability_]:=With[{tt=RandomChoice[transcripts->Range[1,transcripts],((Floor[Total[transcripts]*(1-probability)]/._?Negative->0))]},
(transcripts-ReplacePart[transcripts/._?NumberQ->0,Normal[Counts[tt]]])/._?Negative->0
]


(* ::Input::Initialization:: *)
activeGenome[geneList_,genepos_,grN_,opt_:0]:=
If[opt==0,(*
(* zero mean random and 1 means use exact values*)random selection without competition*)ReplacePart[geneList,grN[[genepos,RandomInteger[{1,Length[grN[[genepos]]]}]]]],
(*Print[geneList," - ",grN];*)
ReplacePart[geneList,grN]
]


(* ::Input::Initialization:: *)
(*competitionTranscripts[ploidy_,grN_,genome_,transcriptList_]:=
With[{task=RandomChoice[grN[[#[[1]]]],Count[transcriptList,#[[1]]]]&/@genome},
Select[Flatten[Table[RandomSample[Select[Flatten[Table[Flatten[(If[Length[Flatten[#]]==0,{},{{#[[1]],#,genome[[i,2]]}}])&/@task[[i]],1],{i,Length[task]}],1],#[[1]]==j[[1]]&]/.{}->{{}},1*ploidy],{j,genome}],1],Length[#]!=0&]]*)


(* ::Input::Initialization:: *)
(*competitionTranscripts[ploidy_,grN_,transcriptList_]:=
With[{task=MapThread[If[#1=={},{},RandomChoice[#1,#2]]&,{grN,transcriptList}]},
(*given a grn eg {{1\[Rule]1,2\[Rule]1,3\[Rule]1},{},{2\[Rule]0,3\[Rule]1}} and transcript list {3,1,1}
we select tasks given the number of transcripts for each gene results in {{2\[Rule]1,1\[Rule]1,1\[Rule]1},{},{2\[Rule]0}} respectively
then we group them by gene
and select based on gene ploidy number
for 1 = {{2\[Rule]1},{1\[Rule]1}}, for 2 ={{2\[Rule]1,2\[Rule]0},{1\[Rule]1,1\[Rule]1}}
*)
RandomSample[#,UpTo[ploidy]]&/@Values[GroupBy[Flatten[task],Keys]]
]*)


(* ::Input::Initialization:: *)
competitionTranscripts[CNV_,grN_,transcriptList_(*,CNV_*)]:=
With[{task=MapThread[If[Flatten[#1]=={}||#2<1,{},RandomChoice[CNV[[Keys[#1]]]->#1,#2](*#1[[RandomChoice[Range[1,Length[#1]]->CNV[[Keys[#1]]],#2]]]*)]&,{grN,transcriptList}]},
(*given a grn eg {{1\[Rule]1,2\[Rule]1,3\[Rule]1},{},{2\[Rule]0,3\[Rule]1}} and transcript list {3,1,1}
we select tasks given the number of transcripts for each gene results in {{2\[Rule]1,1\[Rule]1,1\[Rule]1},{},{2\[Rule]0}} respectively
then we group them by gene
and select based on gene ploidy number
for 1 = {{2\[Rule]1},{1\[Rule]1}}, for 2 ={{2\[Rule]1,2\[Rule]0},{1\[Rule]1,1\[Rule]1}}

*)
(*Print["grn\[Rule] ",grN];
Print["cnv\[Rule] ",CNV];
Print["tl\[Rule] ",transcriptList];
Print["task0\[Rule] ",task];*)


task1=GroupBy[Flatten[task],Keys];
(*Print[Keys[task1],CNV[[Keys[task1]]]];*)
MapThread[RandomSample[#1,UpTo[(*ploidy+*)#2]]&,{Values[task1],CNV[[Keys[task1]]]}]
]


(* ::Input::Initialization:: *)
copiesT[ploidy_,cellactivegeneList_,rnaPolRate_,transcriptionTime_,genome_]:=ploidy*(
Floor[(Total[#*{If[rnaPolRate[[2]]>0,0,1] (*>0 is no cell cycle or gene length*),RandomReal[Sort[{0,rnaPolRate[[1]]}]]}]+rnaPolRate[[2]])&/@(transcriptionTime/If[rnaPolRate[[2]]>0,1&/@genome[[All,2]],genome[[All,2]]])]) *cellactivegeneList(*Active*)


(* ::Input::Initialization:: *)
(*copiesT[cellactivegeneList_,rnaPolRate_,transcriptionTime_,genomeLength_]:=
(*only complete transcripts,rnapol[[2]] is reinitiation distance, and rnapol[[1]] is the rate which is usually 1*)(Floor[((Total[Map[Function[dum,Floor[((transcriptionTime*If[rnaPolRate[[1]]<0,Abs[RandomVariate[NormalDistribution[Abs[rnaPolRate[[1]]],Abs[rnaPolRate[[1]]]]]](*randomly selects a positive rate  with mean=sigma*),1(*rnaPolRate[[1]]*)])-dum)/#]/. _?Negative:> 0 ],If[rnaPolRate[[2]]<=0||rnaPolRate[[2]]<=0.||rnaPolRate[[2]]>#,{0},DeleteCases[Range[0,#,rnaPolRate[[2]]],N[#]]]]])&/@genomeLength)
(**If[rnaPolRate[[1]]>0,Abs[RandomVariate[NormalDistribution[rnaPolRate[[1]],rnaPolRate[[1]]],Length[genomeLength]]]*)]/. _?Negative:> 0)*cellactivegeneList(*Active*)*)


(* ::Input::Initialization:: *)
copiesT[cellactivegeneList_,rnaPolRate_,transcriptionTime_,genomeLength_]:=
(*only complete transcripts,rnapol[[2]] is reinitiation distance, and rnapol[[1]] is the rate which is usually fixed if negative or radom around gausian (mean, signa) *)(Floor[((Total[Map[Function[dum,ratePol=If[rnaPolRate[[1]]>0,(*random rate*)Abs[RandomVariate[NormalDistribution[Abs[rnaPolRate[[1]]],Abs[rnaPolRate[[1]]]]]](*randomly selects a positive rate  with mean=sigma*),(*fixed rate*)Abs[rnaPolRate[[1]]]];Floor[((transcriptionTime*ratePol)-(dum/ratePol))/#]/. _?Negative:> 0 ],If[rnaPolRate[[2]]<=0||rnaPolRate[[2]]<=0.||rnaPolRate[[2]]>#,{0},DeleteCases[Range[0,#,rnaPolRate[[2]]],N[#]]]]])&/@genomeLength)
(**If[rnaPolRate[[1]]>0,Abs[RandomVariate[NormalDistribution[rnaPolRate[[1]],rnaPolRate[[1]]],Length[genomeLength]]]*)]/. _?Negative:> 0)*cellactivegeneList(*Active*)


(* ::Input::Initialization:: *)
transcriptsT[copies_,activegenome_]:=Flatten[Tuples[{#[[1]]},#[[2]]][[1]]&/@Transpose[{activegenome,copies}]]


(* ::Input::Initialization:: *)
randDistofT[transcriptLength_]:=({transcriptLength,0}-RandomInteger[{1,transcriptLength}])


(* ::Input::Initialization:: *)
randDistofT2[transcriptValues_]:=With[{(tt=RandomInteger/@transcriptValues)},{transcriptValues-tt,tt}]


(* ::Input::Initialization:: *)
getCellPhase[cell_,shortestToMphase_]:=Module[{leftcc,a},
(*if [[4,3]] < [[4,4]] then that means the cell has passed enough time to divide so it enters Mphase
if [[4,3]] > [[4,4]] then that means the cell cannot reach Mphase and it will be in a particular cell phase either in G1, S and G2
*)
(*Calculated the amount of time left in the cell cycle*)
leftcc=Subtract@@cell[[4,3;;4]];(*example if [[4,3]] is 2hrs and [[4,4]] is 1.5, [[4,3]]-[[4,4]]=0.5. it means there is only half hour for this cell to divde.*)
If[leftcc>shortestToMphase,a=Select[Drop[Transpose[{Range[1,4],Accumulate[cell[[5,1]]/.{0.->0.,_?Negative->0}]}],-1],#[[2]]<=shortestToMphase&];If[a=={},1,Last[a][[2]]],4]
]


(* ::Input::Initialization:: *)
oneCellTranscriptCompetition[oneCell_,grnmotif_,ploidy_]:=Module[{task,newCell,pp},
(*
grnmotif=transcriptionpara[[1]];
ploidy=transcriptionPara[[3]];
*)
newCell=oneCell;
(*Print["pl->",ploidy];*)
If[(Flatten[newCell[[7,1]]]!={}||Union[newCell[[7,1]]]!={0}),(*transcript list is newtissue[[icell,7,1]]*)(*for each copy of the transcript we pick a task*)
task=cellRelatedFunctions`competitionTranscripts[ploidy,grnmotif,newCell[[7,1]](*,Total[newCell[[6,2]]]*)];
(*the tasks chosen also reflect ploidy level*)
(*Print["task1\[Rule] ",task];*)
pp=.;
(*if there are not enough transcripts to cover all plodiy level we fill the lists
eg for ploidy two we obtained {{3\[Rule]1,3\[Rule]1},{2\[Rule]0}} and that will cause  a conflict if we transpose so after padding we get {{3\[Rule]1,3\[Rule]1},{2\[Rule]0,{}}}
*)
If[task!={},
task=Transpose[PadRight[#,Max[ploidy],pp]&/@task/.pp:>{}]; 
(*Print["task2\[Rule] ",task];
Print["nc\[Rule] ",newCell[[6,2]]];*)
(*Print[{newCell[[6,2]],task}];*)
newCell[[6,2]]=MapThread[cellRelatedFunctions`activeGenome[#1,0,Flatten[#2],1(* zero mean random and 1 means use exact values*)]&,{newCell[[6,2]],task}];(*turns on or off genes based transcripts and competition 
and considers ploidy
eg for a n=2 {{1,1,1},{1,1,1}} and tasks {{3\[Rule]1,2\[Rule]0},{3\[Rule]1,{}}}
will result in a new genome {{1,0,1},{1,1,1}}
*)
];
];

newCell
]


(* ::Input::Initialization:: *)
transcriptDegredation[oneCell_,prob_]:=Module[{newCell},
(*prob=transcriptionPara[[2]] probability of transcripts remaining in the system*)(* maintaining trancript probability, trying to capture and average degredation rate,this may be improved by relating it to length*)
newCell=oneCell;
If[((prob==0)||(newCell[[7,1]]=={})||Union[newCell[[7,1]]]=={0}),
newCell[[7,1]]=newCell[[7,1]]/._?NumberQ:>0;,
If[(prob<1),
newCell[[7,1]]=cellRelatedFunctions`parentalTranscriptsAssociation[newCell[[7,1]],prob];(*random transcripts remain upto a probability returns a list with transcript copy numbers for each gene*)]; 
]; (*if prob=1 then it keeps all the transcripts as is*)
newCell
]


(* ::Input::Initialization:: *)
transcriptCopyNumber[oneCell_,shortestToMphase_,rnapol_,genomeLength_]:=Module[{newCell,transcriptionTime},
newCell=oneCell;
(*calcualte current trancription time available for a cell*)transcriptionTime=Min[{newCell[[4,3]],newCell[[4,4]]+shortestToMphase}];
(*identifies the shortest cell cycle duration*)
(*the shortest cell cycle is the difference between the total cell cycle time and the time that has passed, [[4,3]] is the total cell cycle duration, and [[4,4]] is the time that has relatively passed
example if [[4,3]] is 2hrs and [[4,4]] is 1.5, [[4,3]]-[[4,4]]=0.5. it means there is only half hour for this cell to divde. 
We calcualte this value for evey cell and get the shortest in the tissue. 
All cells that have a value \[LessEqual] shortestToMphase will divide
*)
(*calculates the trancripts possible and their copies based on time and ploidy, we can add an inefficiacy term later*)(*genome = Transpose[{genomeName,genomeLength}], we just need the genomeLength*)
(*active genes = newtissue[[icell,6,2]] {{1,1,1},{1,0,1}...Ploidy}*)
(*transcriptionPara5 rpolII rate 0= no change, 2 is up to 2x faster or -2 slower,i didnt try both fast and slow*)(*places the transcripts in the temp sections until division occurs*)
(*Print["in->",newCell[[7]]];
Print["btn->",newCell[[7,2]]];*)
newCell[[7,2]]=(cellRelatedFunctions`copiesT[Total[newCell[[6,2]]](*activegenes*),rnapol,transcriptionTime,genomeLength]/. _?Negative-> 0);
(*Print["acn->",newCell[[7,2]]];*)
(*Values[AssociationThread[genome[[All,1]],copies]]*);(*temporary vs fixed transcripts, they only become fixed right before division*)
newCell
]


(* ::Input::Initialization:: *)
End[]


(* ::Input::Initialization:: *)
Protect[Evaluate[$Context<>"*"]];


(* ::Input::Initialization:: *)
EndPackage[]
