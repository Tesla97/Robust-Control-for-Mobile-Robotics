# Vehicle-Control

Controllo Robusto di un Pendolo Inverso su Carrello.
Implementazione di un controllore H-infinito robusto , L1 robusto e H2 robusto.

Dalla modellazione matematica del sistema tramite definizione di energia cinetica K e potenziale T (L = K - T Lagrangiana) , il passo successivo è stato
l'implementazione di controllori robusti ottimi (H-inf , H2 , L1) con l'obiettivo di stabilizzare robustamente il processo e ridurre l'effetto sulla dinamica
di eventuali disturbi agenti sul sistema , tramite minimizzazione dei rispettivi guadagni indotti.
Confronto delle prestazioni dei controlli così ottenuti con controllori classici , quali , PID , Controllo Non Lineare , Back-Stepping Method.

Implementazione simulazione in Matlab tramite uso di SimMechanics.

Ideato e Realizzato Interamente da Nicola Corea per il corso di Controllo dei Veicoli LM Automazione.

