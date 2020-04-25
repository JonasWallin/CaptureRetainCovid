# Capture retain model för dödstal

En enkel "capture-retain" model för att prediktera hur många dödsfall kommer vara rapporterad x dagar i framtiden givet att dagen paserat. 
I biologi används det med att man fångar in djur, utan att släppa ut dom, flera gånger. Om man fångar in en betydande del av den underliggande populationen kan skatta populationen med hur sekvensen av skattade faller.

## Model
### Enklaste Capture retain

Modellen har antalet rapporterade dödsfall i följande matris:

``` r
Reported[i,j] - antatelet rapporterade döda för dag i
                raporterat för dag j
```
Vi är intresserade av vad `Reported[i,i+x]` är.
Det modellerena med hjälp av modellen 

``` r
dij = Reported[i,j]  - cumsum(Reported[i,i:(j-1)] 
 # antal nya döda rapporterat dag j
nij = Reported[i,j+x] - cumsum(Reported[i,i:(j-1)]
 # antal möjliga som kan rapporteras

dbin(  x = dij, 
	 size = nij,
	 prob = P[j])

```

Nyckeln är alltså att efter `j` förväntas vi ha lika stor sannolikhet att observa varje möjligt fall `nij`.

Denna modell skattas i filen *MultinomialSimple.R*.


### logit Capture retain
Det är lätt att se att antalet nya fall som rapporterats är mycket färre på helger detta korrigerar vi genom

``` r
P[j] = 1/(1+ exp(- beta[j] + beta.helg * helg))
``` 
Denna modell skattas i filen *MultinomialLogit.R*.



### beta-multinomial capture retain
Ofta kan binomial/multinominal fördelningen ge för säkra skattningar pga att vi inte tar hänsyn okända faktorer. För att korrigera för detta använder istället beta-multionimal fördelningen. Istället för ett fixt sannolikhet, `P[i,j]` låter vi istället `P[i,j]` vara beta fördelat med parameterna `alpha` och `beta`. Vi lägger även till kovariater igenom 

``` r
alpha_j = exp(- alpha[j] + alpha.helg * helg)
beta_j  = exp(- beta[j]   + beta.helg * helg)
#P_j | alpha,beta    ~  dbeta(P[j], alpha_j, beta_j) 
# (dij, nij)| P_j    ~  dbin( dij, nij,  P_j)
``` 

Denna modell skattas i filen *MultinomialBeta.R*.



## Result
Om vi skattar antalet döda efter 14 dagar får vi följande resultat

| model        | loglikelhood           | parameterar  |
| ------------- |:-------------:| -----:|
| beta-Multinomial  (helg)    | -278.4706   |30 |
| beta-Multinomial      | -300.2816   |26 |
| Multinomial  (helg) | -359.4901 | 15|
| Multinomial   | -443.9738 | 13|


## Prediction
För att göra prediktion så använder jag profile liklihood på antalelt döda (uniform prior om Bayesian). 
Man kan också se i prediktionerna varför helg att ta hänsyn till helg är mycket viktigt:
Utan helg effekt:
![Image description](https://github.com/JonasWallin/CaptureRetainCovid/raw/master/data/dag_14_2020-04-24bM_helgfri.jpeg)
Med helg effekt:
![Image description](https://github.com/JonasWallin/CaptureRetainCovid/raw/master/data/dag_14_2020-04-24bM.jpeg)

### data och kod
Datan är tagen från folkhälsomyndigeten, och koden för att samla in den plus skapande av figurer är tagen från [https://github.com/adamaltmejd/covid](https://github.com/adamaltmejd/covid)
