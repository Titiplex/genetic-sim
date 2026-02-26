# Projet : EvoLife 3D — Écosystème évolutif à partir de génomes “sans sémantique”

## Idée générale

Le projet **EvoLife 3D** vise à construire une simulation 3D temps réel où une population d’organismes évolue, se
développe, se reproduit et interagit dans un environnement.  
La particularité est que les **génomes** des organismes sont des **chaînes alphanumériques** (allèles) qui n’ont *
*aucune signification codée à la main** (“pas de gène vitesse”, “pas de gène couleur”, etc.).

Au lieu de définir des traits explicitement, le programme impose seulement :

- un **substrat** (organismes composés de cellules et liaisons),
- un **monde** (ressources, contraintes, coûts),
- un **interpréteur génome → paramètres** (hash/seed → nombres),
  puis laisse les comportements et “traits” apparaître comme conséquence de la dynamique et de la sélection.

---

## Objectif principal

Créer un **écosystème évolutif** où, à partir de génomes arbitraires :

- des organismes simples (1 cellule) peuvent **grandir en 3D** dans n’importe quelle direction,
- la **forme** et la **structure** émergent par développement (ajout de cellules, liaisons, organisation),
- la population se stabilise ou diverge en niches selon l’environnement,
- des stratégies apparaissent : exploration, exploitation, évitement des toxines, opportunisme, prédation, etc.

En bref : observer l’émergence de diversité (formes + comportements) **sans définir à la main la sémantique des gènes**.

---

## Contraintes (volontaires)

- **Aucune correspondance “gène → trait” codée en dur** :  
  un gène/allèle n’est jamais interprété comme “force”, “vision”, “taille”, etc.
- L’interprétation passe par :
    - des fonctions de hachage,
    - des seeds,
    - des paramètres numériques,
    - des règles générales (physique simplifiée, coûts, flux).
- Les “traits” sont définis au niveau du **substrat** (ce qui est possible), pas au niveau des gènes.

---

## Substrat (ce qui rend l’évolution possible)

### 1) Organismes

Chaque organisme est un assemblage de :

- **cellules** : particules 3D (position locale, vitesse locale, énergie),
- **liaisons** (“springs”) : liens entre cellules avec longueur au repos, rigidité, amortissement,
- un état global : position du centre de masse, vitesse globale, énergie totale, âge.

Ainsi, la “forme” n’est pas un paramètre abstrait, mais la conséquence directe des positions des cellules + forces
internes.

### 2) Développement

La croissance consiste à :

- choisir une direction de croissance,
- placer une nouvelle cellule à proximité,
- connecter la cellule à des voisines (liaisons),
- payer un coût énergétique.

La direction de croissance est calculée à partir de :

- gradients environnementaux (nutriments, lumière, toxines, biomasse),
- structure actuelle de l’organisme (répartition des cellules),
- bases directionnelles pseudo-aléatoires dérivées du génome,
- poids numériques dérivés du génome.

Résultat : un organisme peut se développer **dans n’importe quelle direction** sans direction “préférée” codée en dur.

---

## Environnement (pressions de sélection)

Pour simplifier au départ, l’environnement contient trois facteurs majeurs :

1) **Nutriments (N)**  
   Source d’énergie classique, distribuée spatialement. Favorise les organismes qui savent se placer et rester dans des
   zones riches.

2) **Lumière (L)**  
   Source d’énergie alternative (photosynthèse abstraite), souvent plus présente dans certaines zones. Favorise des
   stratégies différentes : rester en hauteur, grossir pour capter plus, etc.

3) **Toxines (X)**  
   Zones dangereuses qui imposent une perte d’énergie (ou augmentent les coûts). Favorise évitement, résistance, ou
   comportements opportunistes.

**Bonus V2 : Biomasse (B)**  
Quand un organisme meurt, il laisse une biomasse diffuse que d’autres peuvent digérer.  
Cela introduit une chaîne trophique minimale : charognage → prédation → recyclage.

---

## Évolution (variation + sélection)

### Variation

Chaque reproduction produit un enfant avec un génome muté :

- substitution/insertion/suppression de caractères dans les allèles,
- ajout/suppression de gènes (augmentation ou réduction de complexité),
- tout reste “random” au niveau symbolique.

### Sélection

La sélection est implicite :

- survivre exige de maintenir un bilan énergétique positif,
- reproduction nécessite énergie et temps,
- prédation ou digestion offrent des gains possibles mais risqués.

Donc les génomes qui engendrent des organismes “fonctionnels” dans cet environnement se propagent naturellement.

---

## Pourquoi ce projet est intéressant

- **Exploration de l’émergence** : voir apparaître des comportements et formes sans règles spécifiques par trait.
- **Life-like simulation** : reproduction, mort, compétition, niches, diversité.
- **Pont entre génétique abstraite et morphogenèse** :  
  le génome ne décrit pas “des traits”, il décrit un système de paramètres qui gouverne un processus.
- **Plateforme extensible** : ensuite on peut ajouter symbiose, organites, multicellularité avancée, communication,
  spécialisation cellulaire, etc.

---

## Finalité du projet (livrable attendu)

Un programme (C++ / OpenGL) qui :

- affiche un monde 3D en temps réel,
- simule une population d’organismes multicellulaires,
- montre croissance, reproduction, mort, interactions,
- permet d’observer l’évolution sur le long terme,
- fournit des métriques simples (population, taille moyenne, biomasse, taux de prédation).

L’objectif n’est pas de reproduire fidèlement la biologie, mais de construire un **substrat minimal** où des phénomènes
biologiques plausibles (compétition, stratégies, morphologies) peuvent émerger.

---

## Extensions envisagées (après la V2)

- vraie grille 3D de champs (diffusion numérique),
- collisions inter-organismes au niveau cellulaire,
- différenciation cellulaire via un GRN/CTRNN local,
- symbiose (engulf → hôte/symbionte),
- instanced rendering + optimisation SoA pour très grandes populations.

## Proposition (10 changements) + implémentation

> Les 10 changements ci-dessous sont **proposés et implémentés dans cette itération** pour rapprocher la simulation d'une dynamique plus plausible biologiquement, moins déterministe, et plus stable à long terme.

### Changements très importants (impact majeur)

1. **Ajout d'un champ d'oxygène environnemental (O2) + gradient** : l'oxygène varie selon altitude, immersion, humidité, biomasse et climat.
2. **Ajout d'un champ de salinité + coût osmotique** : introduit un stress physiologique réaliste et des compromis écologiques.
3. **Plasticité phénotypique asynchrone** : le phénotype n'est plus recalculé à chaque tick, mais selon une cadence adaptative (stress/apprentissage), ce qui améliore réalisme + performance.
4. **Régulation démographique par capacité de charge dynamique** : la population cible dépend désormais du climat/saisons au lieu d'être fixe.

### Autres changements importants

5. **Saisonnalité de la reproduction** : une fenêtre reproductive fluctuante évite une reproduction uniforme non réaliste.
6. **Réduction de l'endogamie** : pénalité de choix du partenaire si compatibilité génétique trop élevée.
7. **Investissement parental héritable** : réserve énergétique parentale transférée au nouveau-né + coût de cooldown modulé.
8. **Dispersion adaptative** : pression de dispersion augmentée par densité locale et trait génétique dédié.
9. **Recyclage biogéochimique** : dépôt local de nutriments (via excrétion/activité) qui rétroagit sur les niches.
10. **Métriques enrichies** : stress oxygène, stress salinité, taux reproductif, complexité génomique moyenne exportés dans les métriques.

### Résultat attendu

- Plus de **contraintes multi-factorielles** (énergie + maladie + oxygène + osmorégulation + climat).
- Des **niches écologiques** plus nettes (zones humides, zones salines, hauteurs oxygénées, etc.).
- Une dynamique évolutive moins “plate” grâce à la saisonnalité, la capacité de charge dynamique et l'investissement parental.
- Un coût CPU mieux maîtrisé quand la population grandit (plasticité phénotypique asynchrone).
