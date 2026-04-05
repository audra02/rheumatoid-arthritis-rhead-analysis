# Rheumatoid arthritis rhead analysis

Šiame repozitorijoje pateikta mūsų grupės DNR metilinimo duomenų analizė, atlikta naudojant `rhead` duomenų rinkinį. Darbo tikslas buvo aprašyti duomenų rinkinį, įvertinti jo kokybę, ieškoti išskirčių, atlikti klasterizavimą, heatmap analizę ir apžvalginę analizę.

## Naudoti duomenys

Analizei pasirinktas duomenų rinkinys iš publikacijos:

**Brooke, O. M., et al.**  
*Rheumatoid arthritis naive T cells share hypermethylation sites with synoviocytes.*  
Arthritis & Rheumatology, 69(3), 550–559, 2017.

Duomenyse pateikiami DNR metilinimo profiliai iš periferinio kraujo mėginių, suskirstytų pagal skirtingus imuninių ląstelių tipus. Tyrime lyginami reumatoidiniu artritu sergantys ir kontroliniai mėginiai.

## Projekto struktūra

Repozitorijoje naudojama tokia struktūra:

- `data/` – duomenų failai
- `plots/` – sugeneruoti grafikai
- `reports/` – lentelės ir kiti ataskaitos failai
- `scripts/` – analizės skriptai
- `README.md` – projekto aprašymas
- `rheumatoid-arthritis-rhead-analysis.Rproj` – RStudio projekto failas

## Komandos nariai ir jų dalys

Projektą atliko:
- **Audra** – 4 ir 5 užduotys
- **Daniel** – 3 ir 8 užduotys
- **Justinas** – 6, 7 ir 9 užduotys


## Skriptai

### `scripts/Daniel.R`
Šiame skripte atlikta:
- duomenų rinkinio aprašomoji analizė
- pagrindiniai aprašomieji grafikai
- heatmap sudarymas iš labiausiai kintančių CpG pozicijų


### `scripts/4_5_QC.R`
Šiame skripte atlikta kokybės kontrolės analizė:
- metilinimo pasiskirstymo įvertinimas skirtinguose CpG regionuose
- mėginių panašumo analizė pagal ląstelių tipus
- PCA analizė


### `scripts/Justinas_(6+7).R`
Šiame skripte atlikta:
- išskirčių paieška taikant Inter-Array Correlation (IAC)
- filtruoto duomenų rinkinio sudarymas
- hierarchinis klasterizavimas
- jautrumo analizė pašalinus vieną mėginį

### `scripts/Justinas_9.R`
Šiame skripte atlikta papildoma apžvalginė analizė:
- mėginių pasiskirstymo pagal plokšteles analizė
- galimo batch effect įvertinimas
- techninių veiksnių aptarimas


## Trumpas projekto rezultatų apibendrinimas

Analizė parodė, kad pagrindinis DNR metilinimo duomenų variacijos šaltinis yra **ląstelių tipas**, o ne diagnozė. Tai matyti tiek iš koreliacijos analizės, tiek iš PCA, klasterizavimo ir heatmap rezultatų.

Išskirčių paieškos metu buvo identifikuoti ir pašalinti techninės kilmės mėginiai, todėl dalis vėlesnių analizių buvo atliktos naudojant filtruotą duomenų rinkinį. Heatmap analizė parodė aiškų mėginių grupavimąsi pagal ląstelių tipą. Papildoma apžvalginė analizė leido pastebėti, kad kai kuriose plokštelėse mėginiai pagal diagnozę pasiskirstę netolygiai, todėl galimas batch effect turi būti vertinamas atsargiai interpretuojant rezultatus.

## GitHub repozitorija

Repozitorija:  
<https://github.com/audra02/rheumatoid-arthritis-rhead-analysis/tree/main>
