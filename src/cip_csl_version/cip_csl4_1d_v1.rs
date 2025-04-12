/// Applique une advection conservatrice CIP–CSL4 en 1D sur le tableau de densités.
///
/// - `density`: tableau mutable des moyennes de densité par cellule.
/// - `u`: vitesse (constante, positive).
/// - `dt`: pas de temps.
/// - `dx`: taille de cellule (nous supposons dx = 1.0 dans ce code).
///
/// La méthode reconstruit, pour chaque cellule i (avec un stencil complet),
/// un polynôme P(x) de degré 4 tel que l'intégrale sur les 5 cellules du stencil
/// reproduit les moyennes données. Puis, on met à jour la densité de la cellule i
/// en intégrant P(x) sur l'intervalle décalé [–0.5+α, 0.5+α] avec α = u·dt.
pub fn advect_cip_csl4_1d_full(density: &mut [f32], u: f32, dt: f32, dx: f32) {
    let n = density.len();
    let alpha = u * dt / dx; // déplacement fractionnaire (dx=1.0 ici)
    let mut new_density = density.to_vec();

    // Construction de la matrice A du système constant (pour un dx=1, pour k = -2,...,2)
    // Pour chaque k dans {-2,-1,0,1,2} la ligne correspondante est :
    // [ I0(k), I1(k), I2(k), I3(k), I4(k) ]
    // où I_n(k) = ∫_(k-0.5)^(k+0.5) x^n dx.
    let mut A = [[0.0f32; 5]; 5];
    for (row, k) in (-2..=2).enumerate() {
        let a = k as f32 - 0.5;
        let b = k as f32 + 0.5;
        A[row][0] = b - a; // = 1 toujours
        A[row][1] = (b*b - a*a) / 2.0;
        A[row][2] = (b*b*b - a*a*a) / 3.0;
        A[row][3] = (b.powi(4) - a.powi(4)) / 4.0;
        A[row][4] = (b.powi(5) - a.powi(5)) / 5.0;
    }

    // Pour chaque cellule i ayant un stencil complet (i de 2 à n-3)
    for i in 2..(n-2) {
        // Construction du vecteur b_vec des moyennes sur le stencil de 5 cellules
        // b_vec[j] = density[i - 2 + j] pour j=0,...,4
        let mut b_vec = [0.0f32; 5];
        for j in 0..5 {
            b_vec[j] = density[i + j - 2];
        }
        // Copie locale de A dans une matrice M (5x5) qui sera modifiée par élimination gaussienne.
        let mut M = A;
        let mut coeff = [0.0f32; 5]; // Coefficients du polynôme à déterminer.

        // Résolution du système M * coeff = b_vec par élimination gaussienne.
        // Phase d'élimination
        for k in 0..5 {
            let pivot = M[k][k];
            for j in k..5 {
                M[k][j] /= pivot;
            }
            b_vec[k] /= pivot;
            for i2 in (k+1)..5 {
                let factor = M[i2][k];
                for j in k..5 {
                    M[i2][j] -= factor * M[k][j];
                }
                b_vec[i2] -= factor * b_vec[k];
            }
        }
        // Retour-substitution
        for k in (0..5).rev() {
            coeff[k] = b_vec[k];
            for j in (k+1)..5 {
                coeff[k] -= M[k][j] * coeff[j];
            }
        }
        // Maintenant, le polynôme reconstruit est :
        // P(x) = coeff[0] + coeff[1]*x + coeff[2]*x^2 + coeff[3]*x^3 + coeff[4]*x^4,
        // défini pour x ∈ [-0.5, 0.5] dans la cellule i.

        // Pour la mise à jour conservatrice, on calcule l'intégrale de P(x) sur l'intervalle
        // [ -0.5 + alpha, 0.5 + alpha ] (le domaine de départ décalé).
        let antideriv = |x: f32| -> f32 {
            coeff[0] * x +
                coeff[1] * x.powi(2) / 2.0 +
                coeff[2] * x.powi(3) / 3.0 +
                coeff[3] * x.powi(4) / 4.0 +
                coeff[4] * x.powi(5) / 5.0
        };
        let flux = antideriv(0.5 + alpha) - antideriv(-0.5 + alpha);
        new_density[i] = flux;
    }

    // Pour les cellules en bord, on laisse la valeur inchangée
    for i in 0..2 {
        new_density[i] = density[i];
    }
    for i in (n-2)..n {
        new_density[i] = density[i];
    }

    // Mise à jour finale du tableau de densité
    density.copy_from_slice(&new_density);
}



/// Fonction de test pour l'advection CIP–CSL4 en 1D.
pub fn test_advect_cip_csl4_1d() {
    // Création d'un vecteur de densités sur 20 cellules
    let n = 20;
    let mut density = vec![0.0f32; n];

    // Initialisation : on crée un "blob" de densité = 1.0 pour les cellules 8 à 12
    for i in 8..=12 {
        density[i] = 1.0;
    }
    for i in 4..=8 {
        density[i] = 0.5;
    }
    for i in 12..=16 {
        density[i] = 0.5;
    }

    let total_mass_before: f32 = density.iter().sum();
    println!("Total mass before advection: {}", total_mass_before);
    println!("Densities before: {:?}", density);

    let u = 1.0;    // vitesse constante
    let dt = 0.1;   // pas de temps
    let dx = 1.0;   // taille de cellule

    // Appel de l'advection conservatrice CIP–CSL4 en 1D
    advect_cip_csl4_1d_full(&mut density, u, dt, dx);

    let total_mass_after: f32 = density.iter().sum();
    println!("Total mass after advection: {}", total_mass_after);
    println!("Densities after advection: {:?}", density);
}

