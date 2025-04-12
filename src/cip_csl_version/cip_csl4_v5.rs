/// Calcule l'inverse d'une matrice 5x5 par élimination gaussienne.
/// On suppose que la matrice est inversible.
fn invert_matrix_5x5(A: [[f32; 5]; 5]) -> [[f32; 5]; 5] {
    let mut M = A;
    let mut inv = [[0.0f32; 5]; 5];
    // Initialise inv en matrice identité
    for i in 0..5 {
        inv[i][i] = 1.0;
    }
    // Élimination gaussienne
    for i in 0..5 {
        let pivot = M[i][i];
        for j in 0..5 {
            M[i][j] /= pivot;
            inv[i][j] /= pivot;
        }
        for k in 0..5 {
            if k != i {
                let factor = M[k][i];
                for j in 0..5 {
                    M[k][j] -= factor * M[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }
    inv
}

/// Retourn l'inverse constant de la matrice A utilisée pour le stencil 5 points.
/// Pour dx = 1, avec k = -2,...,2, la matrice A est définie par
/// A[row][n] = I_n(k) = ∫_(k-0.5)^(k+0.5) x^n dx, pour n=0..4, et k = -2,-1,0,1,2.
/// Les valeurs numériques ci-dessous sont calculées une fois pour toutes.
fn inv_A_const() -> [[f32; 5]; 5] {
    // On calcule A pour k = -2,-1,0,1,2.
    // Pour dx=1, nous avons:
    // k = -2: a = -2.5, b = -1.5 => I0 = 1, I1 = -2, I2 ≈ 4.08333, I3 = -8.5, I4 ≈ 18.0125
    // k = -1: a = -1.5, b = -0.5 => I0 = 1, I1 = -1, I2 ≈ 1.08333, I3 = -1.25, I4 ≈ 1.5125
    // k =  0: a = -0.5, b =  0.5 => I0 = 1, I1 =  0, I2 ≈ 0.08333, I3 =  0, I4 ≈ 0.0125
    // k =  1: a =  0.5, b =  1.5 => I0 = 1, I1 =  1, I2 ≈ 1.08333, I3 =  1.25, I4 ≈ 1.5125
    // k =  2: a =  1.5, b =  2.5 => I0 = 1, I1 =  2, I2 ≈ 4.08333, I3 =  8.5, I4 ≈ 18.0125
    let A = [
        [1.0,   -2.0,    4.08333, -8.5,    18.0125],
        [1.0,   -1.0,    1.08333, -1.25,   1.5125],
        [1.0,    0.0,    0.08333,  0.0,    0.0125],
        [1.0,    1.0,    1.08333,  1.25,   1.5125],
        [1.0,    2.0,    4.08333,  8.5,    18.0125],
    ];
    invert_matrix_5x5(A)
}

/// Calcule l'interpolation CIP–CSL4 en 1D avec limiteur et utilisation de l'inverse pré-calculé.
/// - `density` : tableau mutable des moyennes de densité par cellule.
/// - `u` : vitesse (constante, positive).
/// - `dt` : pas de temps.
/// - `dx` : taille de cellule (ici supposé = 1.0).
///
/// La méthode reconstruit pour chaque cellule (avec un stencil complet, i de 2 à n-3)
/// un polynôme de degré 4 tel que son intégrale sur l'intervalle décalé [ -0.5+alpha, 0.5+alpha ]
/// donne la nouvelle moyenne. Un simple limiteur ramène les valeurs négatives à 0 et plafonne
/// les overshoots en se basant sur la moyenne locale.
pub fn advect_cip_csl4_1d_limited(density: &mut [f32], u: f32, dt: f32, dx: f32) {
    let n = density.len();
    let alpha = u * dt / dx; // déplacement fractionnaire
    let inv_A = inv_A_const(); // Inverse pré-calculé de la matrice A (5x5)
    let mut new_density = density.to_vec();

    // Pour chaque cellule avec un stencil complet : i de 2 à n-3
    for i in 2..(n - 2) {
        // Construit le vecteur b (les moyennes sur le stencil)
        let mut b = [0.0f32; 5];
        for j in 0..5 {
            b[j] = density[i + j - 2];
        }
        // Calcul des coefficients du polynôme : coeff = inv_A * b
        let mut coeff = [0.0f32; 5];
        for irow in 0..5 {
            for icol in 0..5 {
                coeff[irow] += inv_A[irow][icol] * b[icol];
            }
        }
        // Définition du polynôme P(x) = coeff[0] + coeff[1]*x + ... + coeff[4]*x^4
        // sur x ∈ [-0.5, 0.5]. Sa primitive est donnée par :
        let antideriv = |x: f32| -> f32 {
            coeff[0] * x
                + coeff[1] * x.powi(2) / 2.0
                + coeff[2] * x.powi(3) / 3.0
                + coeff[3] * x.powi(4) / 4.0
                + coeff[4] * x.powi(5) / 5.0
        };
        // Intégration sur l'intervalle décalé [ -0.5 + alpha, 0.5 + alpha ]
        let flux = antideriv(0.5 + alpha) - antideriv(-0.5 + alpha);

        // Limiteur : contraindre flux à rester dans l'intervalle [min(b), max(b)]
        let local_min = *b.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        let local_max = *b.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        let flux_limited = flux.clamp(local_min, local_max);
        new_density[i] = flux_limited;
    }

    // Pour les cellules en bord, on laisse inchangé
    for i in 0..2 {
        new_density[i] = density[i];
    }
    for i in (n - 2)..n {
        new_density[i] = density[i];
    }

    // Mise à jour finale du tableau de densité
    density.copy_from_slice(&new_density);
}



/// Fonction de test pour l'advection CIP–CSL4 1D avec limiteur.
pub(crate) fn test_advect_cip_csl4_1d_limited() {
    // On crée 20 cellules avec un profil "blob" lissé.
    let n = 20;
    let mut density = vec![0.0f32; n];

    // Exemple de profil lissé : densités variant progressivement.
    // Par exemple, on met un pic lissé de la cellule 8 à 12.
    for i in 8..=12 {
        density[i] = 1.0;
    }
    // Pour avoir un peu de transition, on peut lisser manuellement autour.
    density[7] = 0.5;
    density[13] = 0.5;

    let total_mass_before: f32 = density.iter().sum();
    println!("Total mass before advection: {}", total_mass_before);
    println!("Densities before: {:?}", density);

    let u = 1.0;    // vitesse constante
    let dt = 0.1;   // pas de temps
    let dx = 1.0;   // taille de cellule

    // Appel de l'advection CIP–CSL4 1D avec limiteur
    advect_cip_csl4_1d_limited(&mut density, u, dt, dx);

    let total_mass_after: f32 = density.iter().sum();
    println!("Total mass after advection: {}", total_mass_after);
    println!("Densities after advection: {:?}", density);
}

