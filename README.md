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

Projektą atliko: **Audra Stepanauskaitė**, **Daniel Volčak**, **Justinas Tomkevičius**.

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

### `scripts/task2/00_find_outliers.R`

Šiame skripte atlikta papildoma išskirtinių mėginių analizė.
Iš pradinio `rhead.rds` duomenų rinkinio pirmiausia pašalinami anksčiau nustatyti išskirtiniai mėginiai.
Tada sudaroma dendrograma su mėginių ID, leidžianti vizualiai įvertinti mėginių grupavimąsi ir galimas papildomas išskirtis.
Remiantis dendrograma ir naujais pastebėjimais, atrenkami papildomi galimi išskirtiniai mėginiai.
Po jų pašalinimo sudaroma nauja dendrograma.

### `scripts/task2/01_remove_outliers_and_statistics.R`

Šis skriptas pašalina iš ankstesnių kokybės kontrolės žingsnių nustatytus išskirtinius mėginius.
Kiekvienai CpG pozicijai Wilcoxon testu palyginamos reumatoidinio artrito ir kontrolinės grupės.
Apskaičiuojamos p reikšmės, efekto dydžiai, FDR metodu koreguotos p reikšmės ir reikšmingumo žyma.
Sukuriamas filtruotas annmatrix objektas be išskirtinių mėginių ir rezultatų lentelė tolesnei analizei.
**Šį skriptą reikią paleisti patį pirmą lokaliai laikant originalius rhead.rds duomenis scripts/task2**
**tik tuomet visi sekantys skriptai veiks, jų eiliškumas nesvarbus.**

### `scripts/task2/02_top10_plots.R`

Šiame skripte atrenkamos statistiškai reikšmingos CpG pozicijos, kurių p_adj < 0.05.
Iš jų pasirenkama TOP10 CpG pozicijų pagal didžiausią absoliutų efekto dydį.
Kiekvienai TOP10 CpG pozicijai sukuriamas boxplot grafikas su atskirais mėginių taškais ir grupės vidurkiu.
Grafikai išsaugomi aplanke plots/task2 ir naudojami vizualiai įvertinti metilinimo skirtumus tarp RA ir kontrolinės grupės.

### `scripts/task2/03_pvalue_histogram.R`

Šiame skripte sudaroma visų CpG pozicijų p reikšmių histograma.
p reikšmės suskirstomos į intervalus nuo 0 iki 1 ir vizualizuojamos naudojant ggplot2.
Histograma naudojama įvertinti bendrą statistinių testų rezultatų pasiskirstymą
ir patikrinti, ar stebimas mažų p reikšmių perteklius, rodo galimus skirtumus tarp RA ir kontrolinės grupės.

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
<https://github.com/audra02/rheumatoid-arthritis-rhead-analysis
