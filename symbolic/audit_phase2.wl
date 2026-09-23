(* ::Package:: *)

(* ============================================================ *)
(* Phase II audit entry point                                   *)
(* ============================================================ *)

(* Purpose:
   Run the authoritative Phase II polynomial derivation audit.

   Mathematical statement:
     Phase2AdmissibleDomain[] && DStarLambda[] == 0
       ==> P14Polynomial[] == 0.

   Command-line use:

     wolframscript -script symbolic/audit_phase2.wl

   Fresh-notebook use:

     Get[FileNameJoin[{projectRoot, "symbolic", "audit_phase2.wl"}]]

   Expected behavior:
     - PASS requires symbolic theorem checks and numerical spot checks.
     - UNRESOLVED symbolic checks count as failure.
     - The reverse implication is not claimed.
*)

ClearAll[
  SymbolicPhase2AuditDirectory,
  SymbolicPhase2PrintSummary,
  RunSymbolicPhase2Audit,
  SymbolicPhase2AuditResult,
  SymbolicPhase2AuditCommandLineQ
];

If[! ValueQ[$AuditVerbosity], $AuditVerbosity = "Summary"];

SymbolicPhase2AuditDirectory[] :=
  If[StringQ[$InputFileName] && $InputFileName =!= "",
    DirectoryName[$InputFileName],
    If[ValueQ[projectRoot],
      FileNameJoin[{projectRoot, "symbolic"}],
      Directory[]
    ]
  ];

SymbolicPhase2PrintSummary[result_Association] := Module[
  {
    checks = result["SymbolicChecks"]["ForwardPreservingChecks"]
  },
  If[$AuditVerbosity === "Quiet", Return[Null]];
  Scan[
    Print["[" <> #["Status"] <> "] " <> #["Name"]] &,
    checks
  ];
  Print[
    "Overall Phase II status: ",
    If[TrueQ[result["PASS"]], "PASS", "FAIL"]
  ];
];

RunSymbolicPhase2Audit[] := Module[
  {
    scriptDirectory,
    logDirectory,
    exportDirectory,
    timestamp,
    result,
    logPath,
    resultPath,
    polynomialPath,
    coefficientPath,
    polynomialTeXPath,
    reportPath
  },

  scriptDirectory = SymbolicPhase2AuditDirectory[];

  Get[FileNameJoin[{scriptDirectory, "verification", "Phase1Verification.wl"}]];
  Get[FileNameJoin[{scriptDirectory, "verification", "Phase2PolynomialDerivation.wl"}]];

  logDirectory = FileNameJoin[{scriptDirectory, "logs"}];
  exportDirectory = FileNameJoin[{scriptDirectory, "exports"}];
  If[! DirectoryQ[logDirectory],
    CreateDirectory[logDirectory, CreateIntermediateDirectories -> True]
  ];
  If[! DirectoryQ[exportDirectory],
    CreateDirectory[exportDirectory, CreateIntermediateDirectories -> True]
  ];

  timestamp =
    DateString[{"Year", "Month", "Day", "_", "Hour", "Minute", "Second"}];

  result = RunPhase2PolynomialDerivation[];
  SymbolicPhase2PrintSummary[result];

  logPath = FileNameJoin[{logDirectory, "phase2_audit_" <> timestamp <> ".txt"}];
  resultPath = FileNameJoin[{logDirectory, "phase2_result_" <> timestamp <> ".wl"}];
  polynomialPath =
    FileNameJoin[{exportDirectory, "P14_polynomial_" <> timestamp <> ".wl"}];
  coefficientPath =
    FileNameJoin[{exportDirectory, "P14_coefficients_" <> timestamp <> ".wl"}];
  polynomialTeXPath =
    FileNameJoin[{exportDirectory, "P14_polynomial_" <> timestamp <> ".tex"}];
  reportPath =
    FileNameJoin[{exportDirectory, "phase2_report_" <> timestamp <> ".md"}];

  Export[
    logPath,
    StringRiffle[
      {
        "Phase II symbolic polynomial-derivation audit",
        "Timestamp: " <> DateString[],
        "Theorem: " <> result["Theorem"],
        "PASS: " <> ToString[result["PASS"]],
        "Forward theorem PASS: " <> ToString[result["ForwardTheoremPASS"]],
        "Numerical spot-check PASS: " <> ToString[result["NumericalSpotCheckPASS"]],
        "P14 is polynomial in s: " <> ToString[result["SymbolicChecks"]["P14IsPolynomialInS"]],
        "Degree in s: " <> ToString[result["P14Degree"]],
        "Coefficient count: " <> ToString[result["P14CoefficientCount"]],
        "Leading coefficient: " <> ToString[result["P14LeadingCoefficient"], InputForm],
        "Forward-preserving checks: " <> ToString[result["SymbolicChecks"]["ForwardPreservingChecks"], InputForm],
        "Spurious-root sources: " <> ToString[result["SpuriousRootSources"], InputForm],
        "Possible degree drops: " <> ToString[result["PossibleDegreeDropConditions"], InputForm],
        "Numerical spot checks: " <> ToString[result["NumericalSpotChecks"], InputForm]
      },
      "\n"
    ],
    "Text"
  ];

  Put[result, resultPath];
  Put[result["P14Polynomial"], polynomialPath];
  Put[result["P14CoefficientList"], coefficientPath];

  Export[
    polynomialTeXPath,
    ToString[TeXForm[result["P14Polynomial"]]],
    "Text"
  ];

  Export[
    reportPath,
    StringRiffle[
      {
        "# Phase II polynomial derivation audit",
        "",
        "The established theorem is one-way:",
        "",
        "$$D^*_\\Lambda(s)=0\\Longrightarrow P_{14}(s)=0.$$",
        "",
        "The reverse implication is not claimed.",
        "",
        "- PASS: " <> ToString[result["PASS"]],
        "- Forward theorem PASS: " <> ToString[result["ForwardTheoremPASS"]],
        "- Numerical spot-check PASS: " <> ToString[result["NumericalSpotCheckPASS"]],
        "- Degree in $s$: " <> ToString[result["P14Degree"]],
        "- Coefficient count: " <> ToString[result["P14CoefficientCount"]],
        "- Leading coefficient: `" <> ToString[result["P14LeadingCoefficient"], InputForm] <> "`",
        "- Polynomial normalization factor: `" <> ToString[result["PolynomialNormalizationFactor"], InputForm] <> "`",
        "- Denominator-clearing eta factor retained in admissibility conditions: `" <> ToString[result["DenominatorClearingEtaFactor"], InputForm] <> "`",
        "",
        "## Forward-preserving checks",
        "",
        "```wolfram",
        ToString[result["SymbolicChecks"]["ForwardPreservingChecks"], InputForm],
        "```",
        "",
        "Spurious roots may be introduced by denominator clearing and by the two squaring operations.",
        "",
        "The full polynomial and coefficient list are exported as Wolfram expressions in `symbolic/exports/`."
      },
      "\n"
    ],
    "Text"
  ];

  Join[
    result,
    <|
      "LogPath" -> logPath,
      "ResultPath" -> resultPath,
      "P14PolynomialPath" -> polynomialPath,
      "P14CoefficientPath" -> coefficientPath,
      "P14TeXPath" -> polynomialTeXPath,
      "ReportPath" -> reportPath,
      "PassFailSummary" ->
        If[TrueQ[result["PASS"]],
          "PASS: Phase II polynomial candidate derivation verified.",
          "FAIL: Phase II polynomial candidate derivation did not satisfy all audit checks."
        ]
    |>
  ]
];

SymbolicPhase2AuditResult = RunSymbolicPhase2Audit[];

SymbolicPhase2AuditCommandLineQ := TrueQ[$FrontEnd === Null];

If[SymbolicPhase2AuditCommandLineQ,
  Print[SymbolicPhase2AuditResult["PassFailSummary"]];
  Print["Log: ", SymbolicPhase2AuditResult["LogPath"]];
  If[TrueQ[SymbolicPhase2AuditResult["PASS"]],
    Exit[0],
    Exit[1]
  ],
  SymbolicPhase2AuditResult
]
