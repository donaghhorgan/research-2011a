(* ::Package:: *)

(* ::Title:: *)
(*Network / fusion center functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions for cooperative energy detection.*)
(*Copyright (C) 2012 Donagh Horgan.*)
(*Email: donaghh@rennes.ucc.ie.*)
(**)
(*This program is free software : you can redistribute it and/or modify*)
(*it under the terms of the GNU General Public License as published by*)
(*the Free Software Foundation, either version 3 of the License, or*)
(*(at your option) any later version.*)
(**)
(*This program is distributed in the hope that it will be useful,*)
(*but WITHOUT ANY WARRANTY; without even the implied warranty of*)
(*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See *)
(*COPYING for more details.*)
(**)
(*You should have received a copy of the GNU General Public License*)
(*along with this program. If not, see http://www.gnu.org/licenses.*)


(* ::Subsection::Closed:: *)
(*Version information*)


(* ::Text:: *)
(*09/08/2012*)
(*1.2*)


(* ::Subsection::Closed:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.3: Added support for n-bit decisions.*)
(*Version 1.2: Added basic support for single bit decision double threshold type networks - more work needs to be done on this to extend support for n-bit decision type networks.*)
(*Version 1.11: Added database lookup and caching functionality to the SampleComplexity function.*)
(*Version 1.1: Introduced RelevantOptions function to simplify option management, restructured ProbabilityOfDetection and SampleComplexity functions so that they carry the same options as their children.*)
(*Version 1.0: First working version, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["Network`"];


(* ::Subsection:: *)
(*Decision probabilities*)


(* ::Subsubsection::Closed:: *)
(*Probability of false alarm (general)*)


ProbabilityOfFalseAlarm;


(* ::Subsubsection::Closed:: *)
(*Probability of acquisition (general)*)


ProbabilityOfAcquisition;


(* ::Subsubsection::Closed:: *)
(*Probability of detection (general)*)


ProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Probability of missed detection (general)*)


ProbabilityOfMissedDetection;


(* ::Subsubsection::Closed:: *)
(*Probability of false alarm (1 bit decision only)*)


FusionCenterProbabilityOfFalseAlarm;


(* ::Subsubsection::Closed:: *)
(*Probability of detecion (1 bit decision only)*)


FusionCenterProbabilityOfDetection;


(* ::Subsection::Closed:: *)
(*Optimum threshold*)


\[Lambda]opt;


(* ::Subsection::Closed:: *)
(*Optimum voting rule*)


k;


(* ::Subsection::Closed:: *)
(*Sample complexity*)


SampleComplexity;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


Needs["AWGN`"];
Needs["Rayleigh`"];
Needs["Nakagami`"];
Needs["Rice`"];
Needs["DBLogging`"];


(* ::Subsection:: *)
(*Decision probabilities*)


(* ::Subsubsection::Closed:: *)
(*Probability of false alarm (general)*)


Options[ProbabilityOfFalseAlarm] = {DecisionBits->\[Infinity],CorrelationCoefficient->0,Method->"Approximate"};
ProbabilityOfFalseAlarm::usage="ProbabilityOfFalseAlarm[M, \[Lambda]] calculates the probability of false alarm for a single energy detector.
ProbabilityOfFalseAlarm[M, \[Lambda], n] calculates the probability of false alarm for a cooperative network of energy detectors with infinite precision decision fusion.
ProbabilityOfFalseAlarm[M, \[Lambda], n, k, DecisionBits->Nb] calculates the probability of false alarm for a cooperative network of energy detectors with Nb bit decision fusion.

The following methods may be specified:

Method\[Rule]\"Approximate\"
Method\[Rule]\"Exact\"

By default, Method\[Rule]\""<>ToString[Method/.Options[ProbabilityOfFalseAlarm]]<>"\".

The number of decision bits used in fusion may be specified by the DecisionBits option. By default, DecisionBits\[Rule]"<>ToString[DecisionBits/.Options[ProbabilityOfFalseAlarm]]<>". DecisionBits \[Element] {1,\[Infinity]}.

Additionally, the average correlation between nodes may be specified with the CorrelationCoefficient option. By default, CorrelationCoefficient\[Rule]"<>ToString[CorrelationCoefficient/.Options[ProbabilityOfFalseAlarm]]<>".";
ProbabilityOfFalseAlarm::k="Error: Must specify a voting rule when DecisionBits\[Rule]`1`";
ProbabilityOfFalseAlarm::Nb="Error: The list `1` must be of length `2` when DecisionBits->`3`.";
ProbabilityOfFalseAlarm::\[Lambda]="Error: The list `1` must be of length `2` when n = `2`.";
ProbabilityOfFalseAlarm[M_,\[Lambda]_,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfFalseAlarm][[All,1]]],Options[target][[All,1]]];
	ProbabilityOfFalseAlarm[M,\[Lambda],n,RelevantOptions[ProbabilityOfFalseAlarm]]
]
ProbabilityOfFalseAlarm[M_,\[Lambda]_,n_,OptionsPattern[]]:=Module[{RelevantOptions, k = Null},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfFalseAlarm][[All,1]]],Options[target][[All,1]]];
	ProbabilityOfFalseAlarm[M,\[Lambda],n,k,RelevantOptions[ProbabilityOfFalseAlarm]]
]
ProbabilityOfFalseAlarm[M_,\[Lambda]_,n_,k_,OptionsPattern[]]:=Module[{RelevantOptions, \[Rho] = OptionValue[CorrelationCoefficient], Nb = OptionValue[DecisionBits]},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfFalseAlarm][[All,1]]],Options[target][[All,1]]];
	Which[
		Nb == \[Infinity],
			Which[
				\[Rho] == 0,
					Which[
						!ListQ[\[Lambda]],
							AWGNProbabilityOfFalseAlarm[M,\[Lambda],n,RelevantOptions[AWGNProbabilityOfFalseAlarm]],
						ListQ[\[Lambda]],
							Which[
								Length[\[Lambda]] == 1,
									AWGNProbabilityOfFalseAlarm[M,Last[\[Lambda]],n,RelevantOptions[AWGNProbabilityOfFalseAlarm]],
								Length[\[Lambda]] == 2,
									AWGNProbabilityOfFalseAlarm[M,Last[\[Lambda]],n,RelevantOptions[AWGNProbabilityOfFalseAlarm]],
								True,
									Undefined
							],
						True,
							Undefined
					],
				True,
					Undefined
			],
		Nb >= 1,
			If[TrueQ[k==Null],Message[ProbabilityOfFalseAlarm::k,Nb];Abort[]];
			Module[{probabilities,temp},
				Which[
					!ListQ[\[Lambda]],
						If[Nb != 1, Message[ProbabilityOfFalseAlarm::Nb, \[Lambda], 2^Nb - 1, Nb]; Abort[]];
						probabilities = AWGNProbabilityOfFalseAlarm[M,\[Lambda],RelevantOptions[AWGNProbabilityOfFalseAlarm]],
					ListQ[\[Lambda]],
						Which[
							!ListQ[First[\[Lambda]]],
								If[Length[\[Lambda]] != 2^Nb - 1, Message[ProbabilityOfFalseAlarm::Nb, \[Lambda], 2^Nb - 1, Nb]; Abort[]];
								temp = Table[AWGNProbabilityOfFalseAlarm[M,x,RelevantOptions[AWGNProbabilityOfFalseAlarm]],{x,\[Lambda]}];
								probabilities = Flatten[{-Differences[Flatten[{1,temp}]],Last[temp]}],
							ListQ[First[\[Lambda]]],
								If[Length[\[Lambda]] != n, Message[ProbabilityOfFalseAlarm::\[Lambda], \[Lambda], n]; Abort[]];
								Table[If[Length[\[Lambda][[i]]] != 2^Nb - 1, Message[ProbabilityOfFalseAlarm::Nb, \[Lambda][[i]], 2^Nb - 1, Nb]; Abort[]],{i,Length[\[Lambda]]}];
								probabilities = Table[temp = Table[AWGNProbabilityOfFalseAlarm[M,x,RelevantOptions[AWGNProbabilityOfFalseAlarm]],{x,\[Lambda][[n0]]}];Flatten[{-Differences[Flatten[{1,temp}]],Last[temp]}],{n0,n}],
							True,
								Undefined
						],
					True,
						Undefined
				];
				FusionCenterProbabilityOfFalseAlarm[probabilities,n,k,\[Rho]]
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Probability of acquisition (general)*)


Options[ProbabilityOfAcquisition] = Options[ProbabilityOfFalseAlarm];
ProbabilityOfAcquisition::usage=StringReplace[ProbabilityOfFalseAlarm::usage,"ProbabilityOfFalseAlarm"->"ProbabilityOfAcquisition"];
ProbabilityOfAcquisition[M_,\[Lambda]_,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfAcquisition][[All,1]]],Options[target][[All,1]]];
	ProbabilityOfAcquisition[M,\[Lambda],n,RelevantOptions[ProbabilityOfAcquisition]]
]
ProbabilityOfAcquisition[M_,\[Lambda]_,n_,OptionsPattern[]]:=Module[{RelevantOptions, k = Null},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfAcquisition][[All,1]]],Options[target][[All,1]]];
	ProbabilityOfAcquisition[M,\[Lambda],n,k,RelevantOptions[ProbabilityOfAcquisition]]
]
ProbabilityOfAcquisition[M_,\[Lambda]_,n_,k_,OptionsPattern[]]:=Module[{RelevantOptions},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfAcquisition][[All,1]]],Options[target][[All,1]]];
	Which[
		!ListQ[\[Lambda]],
			1 - ProbabilityOfFalseAlarm[M,\[Lambda],n,k,RelevantOptions[ProbabilityOfFalseAlarm]],
		ListQ[\[Lambda]],
			Which[
				Length[\[Lambda]] == 1,
					1 - ProbabilityOfFalseAlarm[M,First[\[Lambda]],n,k,RelevantOptions[ProbabilityOfFalseAlarm]],
				Length[\[Lambda]] == 2,
					1 - ProbabilityOfFalseAlarm[M,First[\[Lambda]],n,k,RelevantOptions[ProbabilityOfFalseAlarm]],
				True,
					Undefined
			]
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Probability of detection (general)*)


Options[ProbabilityOfDetection] = {ChannelType->"AWGN",DecisionBits->\[Infinity],CorrelationCoefficient->0,Method->"Approximate",Algorithm->"NGaussian",LargeMN->10,LowSNR->True,Timed->False,MaxTime->600,MaxIterations->1000,DatabaseLookup->False,DatabaseCaching->False};
ProbabilityOfDetection::usage="ProbabilityOfDetection[M, \[Gamma], \[Lambda]] calculates the probability of detection for a single energy detector for the specified channel type.
ProbabilityOfDetection[M, \[Gamma], \[Lambda], n] calculates the probability of detection for a cooperative network of energy detectors with infinite precision decision fusion for the specified channel type.
ProbabilityOfDetection[M, \[Gamma], \[Lambda], n, k, DecisionBits->Nb] calculates the probability of detection for a cooperative network of energy detectors with Nb bit decision fusion for the specified channel type.

The following channel types may be specified:

ChannelType\[Rule]\"AWGN\"
ChannelType\[Rule]\"Rayleigh\"
ChannelType\[Rule]{\"Nakagami\",m}
ChannelType\[Rule]{\"Rice\",K}

By default, ChannelType\[Rule]\""<>ToString[ChannelType/.Options[ProbabilityOfDetection]]<>"\".

For each channel type, the following methods may be specified:

Method\[Rule]\"Approximate\"
Method\[Rule]{\"Approximate\", Algorithm\[Rule]...}
Method\[Rule]\"Exact\"
Method\[Rule]{\"Exact\", Algorithm\[Rule]...}

By default, Method\[Rule]\""<>ToString[Method/.Options[ProbabilityOfDetection]]<>"\". Check the individual channel packages for more information on algorithms.

The number of decision bits transmitted to the fusion center may be specified by the DecisionBits option. By default, DecisionBits\[Rule]"<>ToString[DecisionBits/.Options[ProbabilityOfFalseAlarm]]<>".

Additionally, the average correlation between nodes may be specified with the CorrelationCoefficient option. By default, CorrelationCoefficient\[Rule]"<>ToString[CorrelationCoefficient/.Options[ProbabilityOfDetection]]<>"

Finally, timing and database lookup/caching options may be (exclusively) specified. The timing option is specified by:

Timed\[Rule]"<>ToString[Timed/.Options[ProbabilityOfDetection]]<>"
MaxIterations\[Rule]"<>ToString[MaxIterations/.Options[ProbabilityOfDetection]]<>"
MaxTime\[Rule]"<>ToString[MaxTime/.Options[ProbabilityOfDetection]]<>"

where the above options are the defaults, if not specified. If Timed\[Rule]True, then a {Pd, time} list of values will be returned.

If timing is not used, then database lookup/caching may be enabled. This requires that the DBLogging package be loaded. By default:

DatabaseLookup\[Rule]"<>ToString[DatabaseLookup/.Options[ProbabilityOfDetection]]<>"
DatabaseCaching\[Rule]"<>ToString[DatabaseCaching/.Options[ProbabilityOfDetection]]<>"

and the data is stored in the database specified in DBLogging. For database caching, both the DatabaseLookup and DatabaseCaching options must be set to True.";
ProbabilityOfDetection::opt="`1` and `2` options are mutually exclusive. Aborting...";
ProbabilityOfDetection::k=ProbabilityOfFalseAlarm::k;
ProbabilityOfDetection::Nb=ProbabilityOfFalseAlarm::Nb;
ProbabilityOfDetection::\[Lambda]=ProbabilityOfFalseAlarm::\[Lambda];
ProbabilityOfDetection::\[Gamma]=ProbabilityOfFalseAlarm::\[Lambda];
ProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	ProbabilityOfDetection[M,\[Gamma],\[Lambda],n,RelevantOptions[ProbabilityOfDetection]]
]
ProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,n_,OptionsPattern[]]:=Module[{k=Null},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	ProbabilityOfDetection[M,\[Gamma],\[Lambda],n,k,RelevantOptions[ProbabilityOfDetection]]
]
ProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,n_,k_,OptionsPattern[]]:=Module[{channelType, m, \[Rho] = OptionValue[CorrelationCoefficient], Nb = OptionValue[DecisionBits], RelevantOptions, f, g, time},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	If[OptionValue[DatabaseLookup] && OptionValue[Timed], Message[ProbabilityOfDetection::opt,"DatabaseLookup","Timed"]; Abort[]];
	If[OptionValue[DatabaseCaching] && OptionValue[Timed], Message[ProbabilityOfDetection::opt,"DatabaseCaching","Timed"]; Abort[]];
	If[ListQ[\[Gamma]] && (Length[\[Gamma]] != n), Message[ProbabilityOfDetection::\[Gamma], \[Gamma], n]; Abort[]];
	If[ListQ[OptionValue[ChannelType]],
		If[Length[OptionValue[ChannelType]]==2,
			{channelType,m} = OptionValue[ChannelType],
			channelType = OptionValue[ChannelType][[1]]
		],
		channelType = OptionValue[ChannelType];
	];
	g[\[Gamma]0_,\[Lambda]0_,n0_]:=Switch[channelType,
		"AWGN",
		AWGNProbabilityOfDetection[M,\[Gamma]0,\[Lambda]0,n0,RelevantOptions[AWGNProbabilityOfDetection]],
		"Rayleigh",
		NRayleighProbabilityOfDetection[M,\[Gamma]0,\[Lambda]0,n0,RelevantOptions[NRayleighProbabilityOfDetection]],
		"Nakagami",
		NNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda]0,m,n0,RelevantOptions[NNakagamiProbabilityOfDetection]],
		"Rice",
		NRiceProbabilityOfDetection[M,\[Gamma]0,\[Lambda]0,m,n0,RelevantOptions[NRiceProbabilityOfDetection]],
		_,
		Undefined
	];
	f := Which[
		Nb == \[Infinity],
			Which[
				\[Rho] == 0,
					Which[
						!ListQ[\[Lambda]],
							Which[
								!ListQ[\[Gamma]],
									g[\[Gamma],\[Lambda],n],
								ListQ[\[Gamma]],
									g[Total[\[Gamma]],\[Lambda],n],
								True,
									Undefined
							],
						ListQ[\[Lambda]],
							Which[
								Length[\[Lambda]] == 1,
									Which[
										!ListQ[\[Gamma]],
											g[\[Gamma],Last[\[Lambda]],n],
										ListQ[\[Gamma]],
											g[Total[\[Gamma]],Last[\[Lambda]],n],
										True,
											Undefined
									],
								Length[\[Lambda]] == 2,
									Which[
										!ListQ[\[Gamma]],
											g[\[Gamma],Last[\[Lambda]],n],
										ListQ[\[Gamma]],
											g[Total[\[Gamma]],Last[\[Lambda]],n],
										True,
											Undefined
									],
								True,
									Undefined
							],
						True,
							Undefined
					],
				True,
					Undefined
			],
		Nb >= 1,
			If[TrueQ[k==Null],Message[ProbabilityOfDetection::k,Nb];Abort[]];
			Module[{probabilities,temp},
				Which[
					!ListQ[\[Lambda]],
						If[Nb != 1, Message[ProbabilityOfDetection::Nb, \[Lambda], 2^Nb - 1, Nb]; Abort[]];
						Which[
							!ListQ[\[Gamma]],
								probabilities = g[\[Gamma],\[Lambda],1],
							ListQ[\[Gamma]],
								If[OptionValue[Timed],
									{temp, time} = Table[g[\[Gamma][[n0]],\[Lambda],1],{n0,n}]//Transpose,
									temp = Table[g[\[Gamma][[n0]],\[Lambda],1],{n0,n}]
								];
								probabilities = Table[{1 - temp[[i]], temp[[i]]},{i,n}],
							True,
								Undefined
						],
					ListQ[\[Lambda]],
						Which[
							!ListQ[First[\[Lambda]]],
								If[Length[\[Lambda]] != 2^Nb - 1, Message[ProbabilityOfDetection::Nb, \[Lambda], 2^Nb - 1, Nb]; Abort[]];
								Which[
									!ListQ[\[Gamma]],
										If[OptionValue[Timed],
											{temp, time} = Table[g[\[Gamma],\[Lambda]0,1],{\[Lambda]0,\[Lambda]}]//Transpose,
											temp = Table[g[\[Gamma],\[Lambda]0,1],{\[Lambda]0,\[Lambda]}]
										];
										probabilities = Flatten[{-Differences[Flatten[{1,temp}]],Last[temp]}],
									ListQ[\[Gamma]],
										probabilities = Table[
											If[OptionValue[Timed],
												{temp, time} = Table[g[\[Gamma][[n0]],\[Lambda]0,1],{\[Lambda]0,\[Lambda]}]//Transpose,
												temp = Table[g[\[Gamma][[n0]],\[Lambda]0,1],{\[Lambda]0,\[Lambda]}]
											];
											Flatten[{-Differences[Flatten[{1,temp}]],Last[temp]}],{n0,n}
										],
									True,
										Undefined
								],
							ListQ[First[\[Lambda]]],
								If[Length[\[Lambda]] != n, Message[ProbabilityOfDetection::\[Lambda], \[Lambda], n]; Abort[]];
								Table[If[Length[\[Lambda][[i]]] != 2^Nb - 1, Message[ProbabilityOfDetection::Nb, \[Lambda][[i]], 2^Nb - 1, Nb]; Abort[]],{i,Length[\[Lambda]]}];
								Which[
									!ListQ[\[Gamma]],
										probabilities = Table[
											If[OptionValue[Timed],
												{temp, time} = Table[g[\[Gamma],\[Lambda]0,1],{\[Lambda]0,\[Lambda][[n0]]}]//Transpose,
												temp = Table[g[\[Gamma],\[Lambda]0,1],{\[Lambda]0,\[Lambda][[n0]]}]
											];
											Flatten[{-Differences[Flatten[{1,temp}]],Last[temp]}],{n0,n}
										],
									ListQ[\[Gamma]],
										probabilities = Table[
											If[OptionValue[Timed],
												{temp, time} = Table[g[\[Gamma][[n0]],\[Lambda]0,1],{\[Lambda]0,\[Lambda][[n0]]}]//Transpose,
												temp = Table[g[\[Gamma][[n0]],\[Lambda]0,1],{\[Lambda]0,\[Lambda][[n0]]}]
											];
											Flatten[{-Differences[Flatten[{1,temp}]],Last[temp]}],{n0,n}
										],
									True,
										Undefined
								],
							True,
								Undefined
						],
					True,
						Undefined
				];
				If[OptionValue[Timed],
					{FusionCenterProbabilityOfDetection[probabilities,n,k,\[Rho]],Total[time]},
					FusionCenterProbabilityOfDetection[probabilities,n,k,\[Rho]]
				]
			],
		True,
			Undefined
	];
	If[OptionValue[DatabaseLookup],
		result = GetProbabilityOfDetection[OptionValue[Algorithm],channelType,M,\[Gamma],ProbabilityOfFalseAlarm[M,\[Lambda],n,RelevantOptions[ProbabilityOfFalseAlarm]]//N,n,m];
		If[TrueQ[result==Null],
			result = f;
			If[OptionValue[DatabaseCaching],
				(* For correct retrieval of results later, attempt to convert Pf to rational form *)
				rationalPf = Round[ProbabilityOfFalseAlarm[M,\[Lambda],n,RelevantOptions[ProbabilityOfFalseAlarm]]*10^6//N]/10^6;
				CacheProbabilityOfDetection[OptionValue[Algorithm],channelType,M,\[Gamma],rationalPf,n,m,result//N,OptionValue[Algorithm]];
			];
		];
		result,
		f
	]
]


(* ::Subsubsection:: *)
(*Probability of missed detection (general)*)


Options[ProbabilityOfMissedDetection] = Options[ProbabilityOfDetection];
ProbabilityOfMissedDetection::usage=StringReplace[ProbabilityOfDetection::usage,"ProbabilityOfDetection"->"ProbabilityOfMissedDetection"];
ProbabilityOfMissedDetection[M_,\[Gamma]_,\[Lambda]_,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfMissedDetection][[All,1]]],Options[target][[All,1]]];
	ProbabilityOfMissedDetection[M,\[Gamma],\[Lambda],n,RelevantOptions[ProbabilityOfMissedDetection]]
]
ProbabilityOfMissedDetection[M_,\[Gamma]_,\[Lambda]_,n_,OptionsPattern[]]:=Module[{k=Null},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfMissedDetection][[All,1]]],Options[target][[All,1]]];
	ProbabilityOfMissedDetection[M,\[Gamma],\[Lambda],n,k,RelevantOptions[ProbabilityOfMissedDetection]]
]
ProbabilityOfMissedDetection[M_,\[Gamma]_,\[Lambda]_,n_,k_,OptionsPattern[]]:=Module[{channelType,m,\[Rho]=OptionValue[CorrelationCoefficient],RelevantOptions,f},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[ProbabilityOfMissedDetection][[All,1]]],Options[target][[All,1]]];
	Which[
		!ListQ[\[Lambda]],
			1 - ProbabilityOfDetection[M,\[Gamma],\[Lambda],n,k,RelevantOptions[ProbabilityOfDetection]],
		ListQ[\[Lambda]],
			Which[
				Length[\[Lambda]] == 1,
					1 - ProbabilityOfDetection[M,\[Gamma],First[\[Lambda]],n,k,RelevantOptions[ProbabilityOfDetection]],
				Length[\[Lambda]] == 2,
					1 - ProbabilityOfDetection[M,\[Gamma],First[\[Lambda]],n,k,RelevantOptions[ProbabilityOfDetection]],
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Probability of false alarm (fusion center)*)


FusionCenterProbabilityOfFalseAlarm::usage="FusionCenterProbabilityOfFalseAlarm[Pf, n, k] calculates the overall probability of false alarm for a cooperative network with 1 bit decision fusion.
FusionCenterProbabilityOfFalseAlarm[Pf, n, k, \[Rho]] calculates the overall probability of false alarm for a cooperative network with 1 bit decision fusion and correlated decisions.";
FusionCenterProbabilityOfFalseAlarm[P_,n_?IntegerQ,k_?IntegerQ,\[Rho]_:0]:=FusionCenterProbabilityOfFalseAlarm[{P},n,k,\[Rho]]
FusionCenterProbabilityOfFalseAlarm[P_?ListQ,n_?IntegerQ,k_?IntegerQ,\[Rho]_:0]:=Which[
	Length[P] == 1,
		Which[
			\[Rho] == 0,
				1 - CDF[BinomialDistribution[n, First[P]], k - 1],
			0 < \[Rho] <= 1,
				Which[
					n == 1,
						First[P],
					n == 2,
						Undefined,
					n > 2,
						Sum[Sum[(-1)^i Binomial[l, i] First[P] Product[(\[Rho] (s + 1 - First[P]) + First[P]) / (1 + s \[Rho]),{s, 0, n - l + i - 2}],{i, 0, l}],{l, 0, k - 1}],
					True,
						Undefined
				],
			True,
				Undefined
		],
	Length[P] > 1,
		Which[
			\[Rho] == 0,
				Which[
					!ListQ[First[P]],
						Which[
							Length[P] == 2,
								FusionCenterProbabilityOfFalseAlarm[Last[P],n,k,\[Rho]],
							Length[P] == 3,
								Total[Table[Binomial[n,l] P[[1]]^(n-l) P[[3]]^l / (1-P[[2]])^n,{l,k,n}]],
							TrueQ[N[Mod[Log[2,Length[P]],1]] == 0],
								Module[{decisions=Table[Unique[],{n}],probabilities,Nb=Log[2,Length[P]],p,precision=10000},
									Total[
										Table[
											probabilities=Table[Subscript[p, decisions[[x]]],{x,n}]/.FindInstance[Flatten[{Total[decisions]==x,Table[0<=decisions[[x]]<=2^Nb-1,{x,n}]}],decisions,Integers,precision];
											Table[Product[probabilities[[x,y]],{y,Length[probabilities[[x]]]}],{x,Length[probabilities]}]//Total,
											{x,k,n (2^Nb-1)}
										]/.Table[Subscript[p, x]->P[[x+1]],{x,0,2^Nb-1}]
									]
								],
							True,
								Undefined
						],
					ListQ[First[P]],
						Which[
							TrueQ[N[Mod[Log[2,Length[First[P]]],1]] == 0],
								Module[{decisions=Table[Unique[],{n}],probabilities,Nb=Log[2,Length[First[P]]],p,precision=10000},
									Total[
										Table[
											probabilities=Table[Subscript[p, x, decisions[[x]]],{x,n}]/.FindInstance[Flatten[{Total[decisions]==x,Table[0<=decisions[[x]]<=2^Nb-1,{x,n}]}],decisions,Integers,precision];
											Table[Product[probabilities[[x,y]],{y,Length[probabilities[[x]]]}],{x,Length[probabilities]}]//Total,
											{x,k,n (2^Nb-1)}
										]/.Flatten[Table[Subscript[p, x, d]->P[[x, d + 1]],{x, n},{d, 0, 2^Nb - 1}]]
									]
								],
							True,
								Undefined
						],
					True,
						Undefined
				],
			True,
				Undefined
		],
	True,
		Undefined
]


(* ::Subsubsection::Closed:: *)
(*Probability of detecion (fusion center)*)


FusionCenterProbabilityOfDetection::usage="FusionCenterProbabilityOfDetection[Pd, n, k] calculates the overall probability of detection for a cooperative network with 1 bit decision fusion.
FusionCenterProbabilityOfDetection[Pd, n, k, \[Rho]] calculates the overall probability of detection for a cooperative network with 1 bit decision fusion and correlated decisions.";
FusionCenterProbabilityOfDetection[P_,n_?NumericQ,k_?NumericQ,\[Rho]_:0]:=FusionCenterProbabilityOfFalseAlarm[P,n,k,\[Rho]]


(* ::Subsection:: *)
(*Optimum threshold*)


OptimumParameters[M_?NumericQ,\[Gamma]_,n_]:=Module[{errors,thresholds,\[Lambda]0=Table[{Unique[\[Lambda]]},{n}],k0},
{errors,thresholds}=Transpose[Table[FindMinimum[{ProbabilityOfFalseAlarm[M,\[Lambda]0,n,k,DecisionBits->1]+(1-ProbabilityOfDetection[M,\[Gamma],\[Lambda]0,n,k,DecisionBits->1]),Apply[And,Table[Which[i==1,0<\[Lambda]0[[i,1]]<\[Lambda]0[[i+1,1]],i==n,\[Lambda]0[[i,1]]<\[Infinity],True,\[Lambda]0[[i,1]]<\[Lambda]0[[i+1,1]]],{i,n}]]},Table[{\[Lambda]0[[i,1]],(M(1+\[Gamma]))/2},{i,n}]],{k,n}]];
k0=First[First[Position[errors,Min[errors]]]];
{k0,thresholds[[k0]]}
]


(* ::Subsection::Closed:: *)
(*Optimum voting rule*)


k::usage="k[Pf, Pd, n] calculates the optimum fusion rule for a cooperative network with 1 bit decision fusion.
k[Pf, Pd, n, \[Rho]] calculates the optimum fusion rule for a cooperative network with 1 bit decision fusion and correlated decisions.
k[{\!\(\*SubscriptBox[\(P\), \(00\)]\), \!\(\*SubscriptBox[\(P\), \(01\)]\), ...}, {\!\(\*SubscriptBox[\(P\), \(10\)]\), \!\(\*SubscriptBox[\(P\), \(11\)]\), ...}, n] calculates the optimum fusion rule for a cooperative network with arbitrary decision quantisation.
k[{\!\(\*SubscriptBox[\(P\), \(00\)]\), \!\(\*SubscriptBox[\(P\), \(01\)]\), ...}, {\!\(\*SubscriptBox[\(P\), \(10\)]\), \!\(\*SubscriptBox[\(P\), \(11\)]\), ...}, n, \[Rho]] calculates the optimum fusion rule for a cooperative network with arbitrary decision quantisation and correlated decisions.";
k::length="Error: the list `1` should be the same length as `2`.";
k[Pf_?NumericQ,Pd_?NumericQ,n_?NumericQ,\[Rho]_:0]:=k[{1-Pf,Pf},{1-Pd,Pd},{n,n},\[Rho]]
k[P0_?ListQ,P1_?ListQ,votes_?ListQ,\[Rho]_:0]:=Module[{v,n,Pa,\[CapitalDelta]0,Pf,Pm,\[CapitalDelta]1,Pd},
	If[Length[P0]!=Length[P1],
		Message[k::length,P0,P1];
		Abort[];
	];
	If[\[Rho] == 0,
		Switch[Length[P0],
			2,
			n = votes[[1]];
			{Pa, Pf} = P0;
			{Pm, Pd} = P1;
			Min[n, \[LeftCeiling](n Log[Pm / Pa]) / Log[(Pf Pm) / (Pd Pa)]\[RightCeiling]],
			3,
			{v,n} = votes;
			{Pa, \[CapitalDelta]0, Pf} = P0;
			{Pm, \[CapitalDelta]1, Pd} = P1;
			Min[v, \[LeftCeiling](v Log[Pa (1 - \[CapitalDelta]1) / (Pm (1 - \[CapitalDelta]0))]) / Log[(Pd Pa) / (Pf Pm)]\[RightCeiling]],
			_,
			Module[{a = Table[FusionCenterProbabilityOfFalseAlarm[Pf,n,k,\[Rho]]+1-FusionCenterProbabilityOfDetection[Pd,n,k,\[Rho]],{k,1,n}]},
				Position[a, Min[a]][[1, 1]]
			]
		],
		Module[{a = Table[FusionCenterProbabilityOfFalseAlarm[Pf,n,k,\[Rho]]+1-FusionCenterProbabilityOfDetection[Pd,n,k,\[Rho]],{k,1,n}]},
			Position[a, Min[a]][[1, 1]]
		]
	]
]


(* ::Subsection::Closed:: *)
(*Sample complexity*)


Options[SampleComplexity]={ChannelType->"AWGN",DecisionBits->\[Infinity],CorrelationCoefficient->0,Method->"Approximate",LowSNR->True,Tolerance->10^-6,DatabaseLookup->False,DatabaseCaching->False};
SampleComplexity::usage="SampleComplexity[\[Gamma], Pf, Pd] calculates the number of samples required for a single energy detector to operate with the specified decision probabilities at the given signal to noise ratio.
SampleComplexity[\[Gamma], Pf, Pd, n] calculates the number of samples required for a cooperative network of energy detectors to operate with the specified decision probabilities at the given signal to noise ratio.

The following channel types may be specified:

ChannelType\[Rule]\"AWGN\"
ChannelType\[Rule]\"Rayleigh\"
ChannelType\[Rule]{\"Nakagami\",m}

By default, ChannelType\[Rule]\""<>ToString[ChannelType/.Options[SampleComplexity]]<>"\".

The following methods may be specified for Rayleigh and Nakagami channels:

Method\[Rule]\"Approximate\"
Method\[Rule]\"Exact\"

By default, Method\[Rule]\""<>ToString[Method/.Options[SampleComplexity]]<>"\".

Additionally, the number of decision bits used in fusion may be specified by the DecisionBits option. By default, DecisionBits\[Rule]"<>ToString[DecisionBits/.Options[SampleComplexity]]<>". DecisionBits \[Element] {1,\[Infinity]}.

If 1 bit decision fusion is specified, then the average correlation between nodes may be specified with the CorrelationCoefficient option. By default, CorrelationCoefficient\[Rule]"<>ToString[CorrelationCoefficient/.Options[SampleComplexity]]<>"

Numerical tolerance may be specified using the Tolerance option. By default, Tolerance\[Rule]"<>ToString[Tolerance/.Options[SampleComplexity]]<>"

Finally, database lookup/caching may be enabled. This requires that the DBLogging package be loaded. By default:

DatabaseLookup\[Rule]"<>ToString[DatabaseLookup/.Options[SampleComplexity]]<>"
DatabaseCaching\[Rule]"<>ToString[DatabaseCaching/.Options[SampleComplexity]]<>"

and the data is stored in the database specified in DBLogging.  For database caching, both the DatabaseLookup and DatabaseCaching options must be set to True.";
SampleComplexity::tol="The difference between the result `1` and the constraint `2` was greater than the specified tolerance `3`.";
SampleComplexity::opt="`1` and `2` options are mutually exclusive. Aborting...";
SampleComplexity[\[Gamma]_,Pf_,Pd_,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[SampleComplexity][[All,1]]],Options[target][[All,1]]];
	SampleComplexity[\[Gamma],Pf,Pd,n,RelevantOptions[SampleComplexity]]
]
SampleComplexity[\[Gamma]_,Pf_,Pd_,n_,OptionsPattern[]]:=Module[{result,channelType,m,\[Rho]=OptionValue[CorrelationCoefficient],x,y,tol=OptionValue[Tolerance],RelevantOptions,f},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[SampleComplexity][[All,1]]],Options[target][[All,1]]];
	(*If[OptionValue[DatabaseLookup]&&OptionValue[Timed],
		Message[SampleComplexity::opt,"DatabaseLookup","Timed"];
		Abort[];
	];
	If[OptionValue[DatabaseCaching]&&OptionValue[Timed],
		Message[SampleComplexity::opt,"DatabaseCaching","Timed"];
		Abort[];
	];*)
	If[ListQ[OptionValue[ChannelType]],
		If[Length[OptionValue[ChannelType]]==2,
			{channelType,m} = OptionValue[ChannelType],
			channelType = OptionValue[ChannelType][[1]]
		],
		channelType = OptionValue[ChannelType];
	];
	Switch[channelType,
		"AWGN",
		m = Null,
		"Rayleigh",
		m = 1;
	];
	f:=If[OptionValue[DecisionBits]==\[Infinity],
		If[\[Rho]==0,
			Switch[channelType,
				"AWGN",
				AWGNSampleComplexity[\[Gamma],Pf,Pd,n],
				"Rayleigh",
				NRayleighSampleComplexity[\[Gamma],Pf,Pd,n,RelevantOptions[NRayleighSampleComplexity]],
				"Nakagami",
				NNakagamiSampleComplexity[\[Gamma],Pf,Pd,m,n,RelevantOptions[NNakagamiSampleComplexity]],
				"Rice",
				NRiceSampleComplexity[\[Gamma],Pf,Pd,m,n,RelevantOptions[NRiceSampleComplexity]],
				_,
				SampleComplexity[\[Gamma],Pf,Pd,n,RelevantOptions[SampleComplexity]]
			],
			(* No solution for correlated infinite precision fusion *)
			Undefined
		],
		If[\[Rho]!=0&&n==2,
			(* No solution for correlated fusion when n = 2 *)
			Undefined,
			(* Temporarily disable checks - we'll do our own after *)
			Off[FindRoot::lstol,NMinimize::cvmit];
			{x,y} = {x,y}/.(NMinimize[{Abs[FusionCenterProbabilityOfDetection[x,n,k[y,x,n,\[Rho]],\[Rho]]-Pd]+Abs[FusionCenterProbabilityOfFalseAlarm[y,n,k[y,x,n,\[Rho]],\[Rho]]-Pf],0<x<1&&0<y<1},{x,y}][[2]]);
			On[FindRoot::lstol,NMinimize::cvmit];
			If[!(Abs[FusionCenterProbabilityOfDetection[x,n,k[y,x,n,\[Rho]],\[Rho]] - Pd] <= tol),Message[SampleComplexity::tol, FusionCenterProbabilityOfDetection[x,n,k[y,x,n,\[Rho]],\[Rho]]//N, Pd//N, tol//N]];
			If[!(Abs[FusionCenterProbabilityOfFalseAlarm[y,n,k[y,x,n,\[Rho]],\[Rho]] - Pf] <= tol),Message[SampleComplexity::tol, FusionCenterProbabilityOfFalseAlarm[y,n,k[y,x,n,\[Rho]],\[Rho]]//N, Pf//N, tol//N]];
			Switch[channelType,
				"AWGN",
				AWGNSampleComplexity[\[Gamma],y,x],
				"Rayleigh",
				NRayleighSampleComplexity[\[Gamma],y,x,RelevantOptions[NRayleighSampleComplexity]],
				"Nakagami",
				NNakagamiSampleComplexity[\[Gamma],y,x,m,RelevantOptions[NNakagamiSampleComplexity]],
				"Rice",
				NRiceSampleComplexity[\[Gamma],Pf,Pd,m,RelevantOptions[NRiceSampleComplexity]],
				_,
				SampleComplexity[\[Gamma],y,x,n,RelevantOptions[SampleComplexity]]
			]
		]
	];
	If[OptionValue[DatabaseLookup],
		result = GetSampleComplexity[{channelType,m},OptionValue[DecisionBits],n,\[Gamma],Pf,1 - Pd];
		If[TrueQ[result==Null],
			result = f;
			If[OptionValue[DatabaseCaching],
				CacheSampleComplexity[{channelType,m},OptionValue[DecisionBits],n,\[Gamma],Pf,1 - Pd,result];
			];
		];
		result,
		f
	]
]


End[];


EndPackage[];
