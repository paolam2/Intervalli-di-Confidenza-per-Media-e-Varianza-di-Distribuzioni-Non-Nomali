# Intervalli-di-Confidenza-per-Media-e-Varianza-di-Distribuzioni-Non-Nomali

### Tesi di Laurea Triennale in Statistica Gestionale
*Università La Sapienza, Roma - Anno Accademico 2023/2024*

Questa tesi presenta un'analisi comparativa tra i metodi classici e i nuovi **intervalli di confidenza** (CI) proposti dal Prof. J. D. Curto per la stima di media, varianza e loro differenze/rapporti. L'obiettivo primario dello studio è valutare la probabilità di copertura dei CI (con $\alpha = 0.05$) e la loro robustezza in presenza di distribuzioni non normali. Il confronto è stato condotto attraverso estese **simulazioni Monte Carlo**. I risultati ottenuti indicano che gli intervalli di Curto sono ottimi stimatori per varianza e rapporto di varianze, e buoni stimatori per la media. La stima per la differenza di medie è risultata non significativa, lasciando validi tutti gli intervalli in esame.

## Indice
- Intervalli di confidenza per la varianza
    - Stimatori della curtosi 
    - Intervalli di confidenza basati sugli stimatori della curtosi
- Intervalli di confidenza per il rapporto tra varianze
- Intervalli di confidenza per la media
- Intervalli di confidenza per la differenza di medie

## Contenuto del rrepository
- **[Tesi.pdf](Tesi.pdf)**: L'elaborato completo di tesi.
- **[mc_simulazioni.R](mc_simulazioni.R)**: Codice per generare i campioni Monte Carlo e implementare gli intervalli di confidenza.
