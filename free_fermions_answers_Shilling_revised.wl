(* ::Package:: *)

(* Physics 584, Computational Methods,  
   Mark Alford, Jan 2020
*)

(* To edit and run this file within Mathematica,
Start Mathematica
  Open -> navigate to this file and open it

Make changes
Click "Run Package" or "Run All Code" in the top right corner of the window
Make more changes ...

When finished, 
  File-> Save
*)


Print["This is Exercise 1: free_fermions_questions.m"];

e$1[p_,m_] =  Module[{},
(* E: fill in single particle dispersion relation here *) 
Sqrt[p^2+m^2]
];

(* E: Make a plot of the dispersion relation of electrons (mass=0.51 MeV)
for momenta from -5 MeV to +5 MeV *)
Plot[ e$1[p,0.51], {p,-5,+5} ]


number$density[pF_,m_] = Module[{}, pF^3/(3*Pi^2) ];  

(* E: write the energy density function: *)
energy$density[pF_,m_] = Module[{p},
 (* using e$1[p,m] *)
 (* using the "Assumptions" option to tell it that m>0 and pF>0  *)
(1/Pi^2)*Integrate[e$1[p,m]*p^2, {p,0,pF}, Assumptions->{pF>0,m>0}]
];

    
(*  E: write a function that gives the physical value of 
  Fermi momentum for a given chemical potential: *)
pFermi[mu_,m_] = Module[{}, 
Sqrt[mu^2-m^2]
  ];

(* E: Test that for the physical value of the Fermi momentum,
   d e / d n  is mu : *)
If[ FullSimplify[D[energy$density[pFermi[mu,m],m], mu]/D[number$density[pFermi[mu,m],m], mu]- mu ,{mu>m, m>0}]== 0 
(* calculate de/dn - mu. Use Simplify[] and tell it that mu > m and m > 0  *) 
,
 Print["d e / d n == mu   OK"]
,(*else*)
 Print["**** d e / d n == mu failed"]
,(*undefined*)
 Print["!!!! d e / d n == mu undefined"]
];
    
(*  EXERCISE: write a function that gives the pressure: *)
pressure[mu_,m_] = Module[{pF}, (* fill this in *) 
pF = pFermi[mu,m];
FullSimplify[mu* number$density[pF,m]- energy$density[pF,m],{mu>m, m>0}]
];
    
(* As a check, we now evaluate the pressure in MeV^4 of a gas of electrons with 
chemical potential 4 MeV: *)
Print["electron pressure at mu=4MeV is: ", pressure[4.0, 0.51], " MeV^4"];


(* EXERCISE: Test that d pressure / d mu == n: *)


If[ FullSimplify[ D[pressure[mu,m],mu]-number$density[pFermi[mu,m],m],{mu>m, m>0}] == 0
,
 Print["d p / d mu == n   OK"]
,(*else*)
 Print["**** d p / d mu == n failed"]
,(*undefined*)
 Print["!!!! d p / d mu == n undefined"]
];


(* Now we will obtain ultrarelativistic (mu>>m) limits of the number density, 
energy density, and pressure, by re-evaluating 
the general expressions using the ultra-relativistic dispersion relation *)

(* EXERCISE: Ultrarelativistic limit of single particle dispersion relation
 (at such high momenta we can neglect the mass) : *)
e$ultra[p_] = Module[{}, Abs[p] ];

(* EXERCISE: As simple check, make a single plot showing both the exact 
dispersion relation for electrons and the ultrarelativistic approximation to it,
for momentum from -3 MeV to +3 MeV *)
Plot[{e$1[p,0.51],e$ultra[p]},{p,-3,3}]

(* EXERCISE: Now use the ultrarelativistic dispersion relation to generate
ultrarelativistic expressions for the pressure etc: *)
pFermi$ultra[mu_] = mu ;
number$density$ultra[pF_] = Module[{}, pF^3/(3*Pi^2) ];
energy$density$ultra[pF_] = Module[{}, (1/Pi^2)*Integrate[e$ultra[p]*p^2, {p,0,pF}, Assumptions->{pF>0}]];
pressure$ultra[mu_] = Module[{pF}, 
pF = pFermi$ultra[mu];
FullSimplify[mu* number$density$ultra[pF]- energy$density$ultra[pF]]
];

(* EXERCISE: Use the Series[ ] function to check that the ultrarelativistic
expression for the pressure agrees with the m->0 limit of pressure[mu,m]: *)

If[ Normal[Series[pressure[mu,m],{m,0,1},Assumptions->{mu>0}]]-pressure$ultra[mu] == 0
,
 Print["Ultrarelativistic limit of pressure OK"]
,(*else*)
 Print["**** Ultrarelativistic limit of pressure failed!"]
,(*undefined*)
 Print["!!!! Ultrarelativistic limit of pressure undefined"]
];


(* Now we will obtain nonrelativistic limits of the number density, energy
density, and pressure, by assuming (mu-m) << m and re-evaluating the general
expressions using a non-relativistic dispersion relation e$NR and
corresponding expression pFermi$NR for the Fermi momentum: *)

(* non-relativistic single-particle energy, including rest mass: *)
e$NR[p_,m_] = Module[{}, m + p^2/(2*m) ]; 

(* EXERCISE: As simple check, make a single plot showing the exact, 
ultrarelativistic, and nonrelativistic dispersion relations for electrons, 
for momentum from -3 MeV to +3 MeV: *)
Plot[ {e$1[p,0.51],e$ultra[p],e$NR[p,0.51]},{p,-3,3}]


(* EXERCISE: Now use the nonrelativistic dispersion relation to generate
nonrelativistic expressions for the pressure etc: *)
pFermi$NR[mu_,m_]= Module[{}, Sqrt[(2 m)(mu-m)]];
number$density$NR[pF_,m_] = Module[{},  pF^3/(3*Pi^2) ]; (* as before *)
energy$density$NR[pF_,m_] = Module[{}, 
(1/Pi^2)*Integrate[e$NR[p,m]*p^2, {p,0,pF}, Assumptions->{pF>0,m>0}]
];
pressure$NR[mu_,m_] = Module[{pF},  
pF = pFermi$NR[mu,m];
FullSimplify[mu* number$density$NR[pF,m]- energy$density$NR[pF,m],{mu>m, m>0}]
];

If[ Normal[Series[ pressure[mu,m]-pressure$NR[mu,m],{mu,m,3} , Assumptions->{mu>0,m>0}]]  == 0
,
 Print["Nonrelativistic limit of pressure OK"]
,(*else*)
 Print["**** Nonrelativistic limit of pressure failed!"]
,(*undefined*)
 Print["!!!! Nonrelativistic limit of pressure undefined!"]
];

(* EXERCISE: Do the same check for the energy density. 
Since it's a function of pF it is natural to expand in terms of pF. *)

If[ Normal[Series[ energy$density[pF,m]-energy$density[pF,m],{pF,0,3} , Assumptions->{m>0}]] == 0
,
 Print["Nonrelativistic limit of energy density OK"]
,(*else*)
 Print["**** Nonrelativistic limit of energy density failed!"]
,(*undefined*)
 Print["!!!! Nonrelativistic limit of energy density undefined!"]
];

(* EXERCISE: Finally, check that  energy$density$NR is 
 the derivative of pressure$NR with respect to mu: *)
 
If[ FullSimplify[D[pressure$NR[mu,m],mu]- number$density$NR[pFermi$NR[mu,m],m],{mu>0,m>0}]== 0
,
 Print["Nonrelativistic d p / d mu == n   OK"]
,(*else*)
 Print["**** Nonrelativistic d p / d mu == n failed"]
,(*undefined*)
 Print["!!!! Nonrelativistic d p / d mu == n undefined!"]
];












