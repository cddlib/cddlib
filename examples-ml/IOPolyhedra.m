(* IOPolyhedra.m: written by Komei Fukuda;  Version-980111  *)

ReadPolyhedraData[filename_String]:=
	Block[{rf,outlist={},row={},m,n,i,val,st=" ",ty},
		rf=OpenRead[filename];
		While[st=!=EndOfFile && st=!="begin",
			st=Read[rf,String];
		];
		{m,n,ty}=Read[rf,{Number,Number,Word}];
		st=Skip[rf,String];
		Print["m=", m,",   n=",n, " type=", ty];
		Do[
		  If[ty!="rational", 
		  	row=Read[rf,Table[Number,{n}]],
		  	row=Read[rf,Table[Word,{n}]];row=ToExpression[row]
		  ];
		  AppendTo[outlist,row],
		  {i,m}
		];
		Close[rf];
		outlist
	];


ReadIncidenceData[filename_String]:=
	Block[{rf,outlist={},row={},m,n,n1,i,val,st=" ",inc,colon},
		rf=OpenRead[filename];
		While[st=!=EndOfFile && st=!="begin",
			st=Read[rf,String];
		];
		{m,n,n1}=Read[rf,{Number,Number,Number}];
		Print["m=", m,",   n=",n];
		Do[
			{inc}=Read[rf,{Number}];
			colon=Read[rf,{Word}];
			row=Read[rf,Table[Number,{Abs[inc]}]];
		  	If[inc>0,
		  		AppendTo[outlist,row],
		  		AppendTo[outlist,Complement[Range[n1],row]]
		  	],
		  	{i,m}
		];
		Close[rf];
		outlist
	];


ReadAdjacencyData[filename_String]:=
	Block[{rf,outlist={},row={},v,num,i,val,st=" ",deg,colon},
		rf=OpenRead[filename];
		While[st=!=EndOfFile && st=!="begin",
			st=Read[rf,String];
		];
		{v}=Read[rf,{Number}];
		Print["v =", v];
		Do[
			{num, deg}=Read[rf,{Number, Number}];
			colon=Read[rf,{Word}];
			If[deg >0, row=Read[rf,Table[Number,{deg}]], row={}];
		  	AppendTo[outlist,row],
		  	{i,v}
		];
		Close[rf];
		outlist
	];


WritePolyhedraIne[a_List,b_List,filename_String,type_String:"real"]:=
	Block[{wf,m,n,i,j,str},
		wf=OpenWrite[filename];
		{m,n}=Dimensions[a];
		Write[wf,OutputForm["H-representation"]];
		Write[wf,OutputForm["begin"]];
		Write[wf,m,OutputForm["  "],n+1,OutputForm["  "],OutputForm[type]];
		Do[
			str="";
			Do[
				If[j==0,
					str=StringJoin[str,"  ",ToString[b[[i]]]],
					str=StringJoin[str,"  ",ToString[-a[[i,j]],FormatType->InputForm]]
				],
		  	{j,0,n}
		  ];
		  Write[wf,OutputForm[str]],
		  {i,m}
		];
		Write[wf,OutputForm["end"]];
		Write[wf,OutputForm["incidence"]];
		Write[wf,OutputForm["input_incidence"]];
		Write[wf,OutputForm["input_adjacency"]];
		Write[wf,OutputForm["adjacency"]];
		Close[wf];
	];

VEtoIne[a_List,b_List,filename_String]:=WritePolyhedraIne[a,b,filename];

WritePolyhedraExt[a_List,filename_String,type_String:"real"]:=
	Block[{wf,m,n,i,j,str},
		wf=OpenWrite[filename];
		{m,n}=Dimensions[a];
		Write[wf,OutputForm["V-representation"]];
		Write[wf,OutputForm["begin"]];
		Write[wf,m,OutputForm["  "],n+1,OutputForm["  "],OutputForm[type]];
		Do[
			str="";
			Do[
				If[j==0,
					str=StringJoin[str,"  ",ToString[1]],
					str=StringJoin[str,"  ",ToString[a[[i,j]],FormatType->InputForm]]
				],
		  	{j,0,n}
		  ];
		  Write[wf,OutputForm[str]],
		  {i,m}
		];
		Write[wf,OutputForm["end"]];
		Write[wf,OutputForm["hull"]];
		Write[wf,OutputForm["incidence"]];
		Write[wf,OutputForm["input_incidence"]];
		Write[wf,OutputForm["input_adjacency"]];
		Write[wf,OutputForm["adjacency"]];
		Close[wf];
	];
