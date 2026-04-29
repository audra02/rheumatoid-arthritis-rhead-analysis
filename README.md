# Rheumatoid arthritis rhead analysis

Šioje repozitorijoje pateikiama mūsų grupės DNR metilinimo duomenų analizė, atlikta naudojant `rhead` duomenų rinkinį.
Darbo tikslas buvo aprašyti duomenų rinkinį, įvertinti jo kokybę, nustatyti galimas išskirtis, atlikti klasterizavimo, PCA, heatmap ir diferencinio metilinimo analizę.
Taip pat buvo paruošti grafikai, rezultatų lentelės ir failai genų ontologijos analizei.

## Naudoti duomenys

Analizei naudotas duomenų rinkinys iš publikacijos:

**Brooke, O. M., et al.**  
*Rheumatoid arthritis naive T cells share hypermethylation sites with synoviocytes.*  
Arthritis & Rheumatology, 69(3), 550–559, 2017.

Duomenų rinkinyje pateikiami DNR metilinimo profiliai iš reumatoidiniu artritu sergančių ir kontrolinės grupės asmenų mėginių.
Mėginiai priklauso skirtingiems imuninių ląstelių tipams.
Pagrindinis analizės objektas buvo CpG pozicijų metilinimo lygiai ir jų skirtumai tarp reumatoidinio artrito bei kontrolinės grupės mėginių.

## Projekto struktūra

Repozitorijoje naudojama tokia struktūra:

- `plots/` – sugeneruoti grafikai.
- `plots/task1/` – pirmos užduoties grafikai.
- `plots/task2/` – antros užduoties grafikai.
- `reports/` – ataskaitos failai.
- `reports/task1_rhead.docx` – pirmos užduoties ataskaita.
- `scripts/` – analizės skriptai.
- `scripts/task1/` – pirmos užduoties analizės skriptai.
- `scripts/task2/` – antros užduoties analizės skriptai.
- `.gitignore` – failas, nurodantis, kurie failai neturi būti sekami su `git`.
- `README.md` – projekto aprašymas.
- `rheumatoid-arthritis-rhead-analysis.Rproj` – RStudio projekto failas.

## Grupės nariai

Projektą atliko:

- **Audra** – 4 ir 5 užduotys.
- **Daniel** – 3 ir 8 užduotys.
- **Justinas** – 6, 7 ir 9 užduotys.

## Pirmos užduoties skriptai

### `scripts/task1/4_5_QC.R`

Šiame skripte atlikta kokybės kontrolės analizė.
Įvertintas metilinimo lygio pasiskirstymas skirtinguose CpG regionuose.
Atlikta mėginių panašumo analizė pagal ląstelių tipus.
Atlikta PCA analizė, leidžianti įvertinti pagrindinius duomenų variacijos šaltinius.

### `scripts/task1/Daniel.R`

Šiame skripte atlikta duomenų rinkinio aprašomoji analizė.
Paruošti pagrindiniai aprašomieji grafikai.
Sudaryta heatmap diagrama iš labiausiai kintančių CpG pozicijų.

### `scripts/task1/Justinas_(6+7).R`

Šiame skripte atlikta išskirčių paieška taikant Inter-Array Correlation (IAC) metodą.
Sudarytas filtruotas duomenų rinkinys be nustatytų išskirčių.
Atliktas hierarchinis klasterizavimas.
Atlikta jautrumo analizė pašalinus vieną mėginį.

### `scripts/task1/Justinas_9.R`

Šiame skripte atlikta papildoma apžvalginė analizė.
Analizuotas mėginių pasiskirstymas pagal plokšteles.
Įvertinta galima techninių veiksnių, tokių kaip plokštelė ar Sentrix ID, įtaka rezultatams.
Rezultatai panaudoti atsargiau interpretuojant galimą batch effect.

## Antros užduoties skriptai

### `scripts/task2/01_remove_outliers_and_statistics.R`

Šiame skripte pašalinami ankstesnės kokybės kontrolės metu nustatyti išskirtiniai mėginiai.
Kiekvienai CpG pozicijai palyginami reumatoidinio artrito ir kontrolinės grupės mėginiai.
Apskaičiuojamos p reikšmės, efekto dydžiai, metilinimo lygio vidurkiai grupėse ir reikšmingumo žymos.
Sukuriamas filtruotas `annmatrix` objektas ir rezultatų lentelė tolesnei analizei.

### `scripts/task2/02_top10_plots.R`

Šiame skripte atrenkamos CpG pozicijos su mažiausiomis p reikšmėmis.
Paruošiami TOP10 CpG pozicijų grafikai.
Grafikai naudojami vizualiai įvertinti metilinimo skirtumus tarp reumatoidinio artrito ir kontrolinės grupės mėginių.

### `scripts/task2/03_pvalue_histogram.R`

Šiame skripte sudaroma p reikšmių histograma.
Ji naudojama įvertinti bendrą statistinių testų rezultatų pasiskirstymą.
Pagal histogramą galima spręsti, ar duomenyse matomas signalas, susijęs su grupių skirtumais.

### `scripts/task2/04_volcano_plot.R`

Šiame skripte sudaromas volcano plot grafikas.
Grafike pavaizduojamas CpG pozicijų efekto dydis ir statistinis reikšmingumas.
Šis grafikas padeda išskirti CpG pozicijas, kurios turi ir didesnį metilinimo skirtumą, ir mažesnę p reikšmę.

### `scripts/task2/05_manhattan_plot.R`

Šiame skripte sudaromas Manhattan plot grafikas.
Grafike CpG pozicijos vaizduojamos pagal jų vietą genome.
Šis grafikas leidžia įvertinti, ar stipresni signalai telkiasi tam tikrose chromosomose ar genomo srityse.

### `scripts/task2/06_prepare_GO.R`

Šiame skripte paruošiami genų sąrašai genų ontologijos analizei.
`foreground` sąrašui naudojami genai, susiję su patikimiausiais CpG skirtumais.
`background` sąrašui naudojami genai, susiję su visomis tirtomis CpG pozicijomis.
Šie failai gali būti naudojami išoriniame GO analizės įrankyje.

### `scripts/task2/07_Apzvalgine.R`

Šiame skripte atlikta papildoma apžvalginė analizė.
Analizė skirta geriau apibendrinti gautus diferencinio metilinimo rezultatus.
Ji padeda įvertinti, kokie bendri dėsningumai matomi po pagrindinių statistinių palyginimų.

## Trumpas projekto rezultatų apibendrinimas

Analizė parodė, kad pagrindinis DNR metilinimo duomenų variacijos šaltinis yra ląstelių tipas.
Šis dėsningumas matomas koreliacijos analizėje, PCA grafikuose, klasterizavimo rezultatuose ir heatmap diagramoje.
Diagnozės efektas buvo silpnesnis nei ląstelių tipo efektas.
Išskirčių paieškos metu buvo nustatyti ir pašalinti keli techninės kilmės mėginiai.
Po filtravimo dalis vėlesnių analizių buvo atliktos naudojant išvalytą duomenų rinkinį.
Papildoma apžvalginė analizė parodė, kad kai kurie techniniai veiksniai gali turėti įtakos rezultatų interpretacijai.
Dėl to galimas batch effect turi būti vertinamas atsargiai.

Antroje užduotyje kiekvienai CpG pozicijai buvo palyginti reumatoidinio artrito ir kontrolinės grupės mėginiai.
Buvo apskaičiuotos p reikšmės ir metilinimo skirtumai tarp grupių.
Rezultatai pavaizduoti naudojant TOP10 CpG pozicijų grafikus, p reikšmių histogramą, volcano plot ir Manhattan plot.
Taip pat buvo paruošti genų sąrašai tolimesnei genų ontologijos analizei.

## GitHub repozitorija

Repozitorija:  
<https://github.com/audra02/rheumatoid-arthritis-rhead-analysis/tree/main>
