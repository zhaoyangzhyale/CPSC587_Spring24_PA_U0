// matrices should be an array of arrays in row-major format
// [
// [m00, m01, ..., m0n],
// [m10, m11, ..., m1n],
// ...
// [mm0, mm1, ..., mmn]
// ]
export function add_matrix_matrix(m1, m2) {
    if (m1.length !== m2.length || m1[0].length !== m2[0].length) {
        throw new Error("Matrices dimensions must be the same.");
    }

    let result = new Array(m1.length);
    for (let i = 0; i < m1.length; i++) {
        result[i] = new Array(m1[i].length);
        for (let j = 0; j < m1[i].length; j++) {
            result[i][j] = m1[i][j] + m2[i][j];
        }
    }
    return result;
}

// matrices should be an array of arrays in row-major format
// [
// [m00, m01, ..., m0n],
// [m10, m11, ..., m1n],
// ...
// [mm0, mm1, ..., mmn]
// ]
export function sub_matrix_matrix(m1, m2) {
    if (m1.length !== m2.length || m1[0].length !== m2[0].length) {
        throw new Error("Matrices dimensions must be the same.");
    }

    let result = new Array(m1.length);
    for (let i = 0; i < m1.length; i++) {
        result[i] = new Array(m1[i].length);
        for (let j = 0; j < m1[i].length; j++) {
            result[i][j] = m1[i][j] - m2[i][j];
        }
    }
    return result;
}

// matrix should be an array of arrays in row-major format
// [
// [m00, m01, ..., m0n],
// [m10, m11, ..., m1n],
// ...
// [mm0, mm1, ..., mmn]
// ]
export function frobenius_norm_matrix(m) {
    let sum = 0;

    for (let i = 0; i < m.length; i++) {
        for (let j = 0; j < m[i].length; j++) {
            sum += m[i][j] * m[i][j];
        }
    }

    // Return the square root of the sum
    return Math.sqrt(sum);
}

// matrix should be an array of arrays in row-major format
// [
// [m00, m01, ..., m0n],
// [m10, m11, ..., m1n],
// ...
// [mm0, mm1, ..., mmn]
// ]
export function mul_matrix_scalar(m, scalar) {
    let result = new Array(m.length);

    for (let i = 0; i < m.length; i++) {
        result[i] = new Array(m[i].length);

        for (let j = 0; j < m[i].length; j++) {
            result[i][j] = m[i][j] * scalar;
        }
    }

    return result;
}

// matrix should be an array of arrays in row-major format
// [
// [m00, m01, ..., m0n],
// [m10, m11, ..., m1n],
// ...
// [mm0, mm1, ..., mmn]
// ]
export function div_matrix_scalar(m, scalar) {
    let result = new Array(m.length);

    for (let i = 0; i < m.length; i++) {
        result[i] = new Array(m[i].length);

        for (let j = 0; j < m[i].length; j++) {
            result[i][j] = m[i][j] / scalar;
        }
    }

    return result;
}

// matrix should be an array of arrays in row-major format
// [
// [m00, m01, ..., m0n],
// [m10, m11, ..., m1n],
// ...
// [mm0, mm1, ..., mmn]
// ]
export function normalized_matrix(m) {
    let f = frobenius_norm_matrix(m);
    return div_matrix_scalar(m, f);
}

// vectors should be an array [v0, v1, v2, ..., vn]
export function add_vector_vector(v1, v2) {
    if (v1.length !== v2.length) {
        throw new Error("Vectors must be of the same length.");
    }

    let result = new Array(v1.length);
    for (let i = 0; i < v1.length; i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

// matrices should be an array of arrays in row-major format
// [
// [m00, m01, ..., m0n],
// [m10, m11, ..., m1n],
// ...
// [mm0, mm1, ..., mmn]
// ]
export function mul_matrix_matrix(m1, m2) {
    if (m1[0].length !== m2.length) {
        throw new Error('Incompatible matrix dimensions');
    }

    const result = new Array(m1.length).fill(0).map(() => new Array(m2[0].length).fill(0));

    for (let i = 0; i < m1.length; i++) {
        for (let j = 0; j < m2[0].length; j++) {
            for (let k = 0; k < m1[0].length; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }

    return result;
}

// matrix should be an array of arrays in row-major format
// [
// [m00, m01],
// [m10, m11]
// ]
// vector should be in format
// [ v0, v1 ]
export function mul_matrix_2x2_vector_2x1(matrix, vector) {
    let res = mul_matrix_matrix( matrix, [[vector[0]], [vector[1]]] );
    return [res[0][0], res[1][0]];
}

// matrix should be an array of arrays in row-major format
// [
// [m00, m01, m02],
// [m10, m11, m12],
// [m20, m21, m22]
// ]
// vector should be in format
// [ v0, v1, v2 ]
export function mul_matrix_3x3_vector_3x1(matrix, vector) {
    let res = mul_matrix_matrix(matrix, [[vector[0]], [vector[1]], [vector[2]]]);
    return [res[0][0], res[1][0], res[2][0]];
}

// matrix should be an array of arrays in row-major format
// [
// [m00, m01, m02],
// [m10, m11, m12],
// [m20, m21, m22]
// ]
// vector should be in format
// [ v0, v1 ]
export function mul_matrix_3x3_vector_2x1(matrix, vector, pad_value_at_end_of_vector=1.0) {
    let res = mul_matrix_matrix(matrix, [[vector[0]], [vector[1]], [pad_value_at_end_of_vector]]);
    return [res[0][0], res[1][0]];
}

// matrix should be an array of arrays in row-major format
// [
// [m00, m01, m02, m03],
// [m10, m11, m12, m13],
// [m20, m21, m22, m23],
// [m30, m31, m32, m33],
// ]
// vector should be in format
// [ v0, v1, v2, v3 ]
export function mul_matrix_4x4_vector_4x1(matrix, vector) {
    let res = mul_matrix_matrix(matrix, [[vector[0]], [vector[1]], [vector[2]], [vector[3]]]);
    return [res[0][0], res[1][0], res[2][0], res[3][0]];
}

// matrix should be an array of arrays in row-major format
// [
// [m00, m01, m02, m03],
// [m10, m11, m12, m13],
// [m20, m21, m22, m23],
// [m30, m31, m32, m33],
// ]
// vector should be in format
// [ v0, v1, v2 ]
export function mul_matrix_4x4_vector_3x1(matrix, vector, pad_value_at_end_of_vector=1.0) {
    let res = mul_matrix_matrix(matrix, [[vector[0]], [vector[1]], [vector[2]], [vector[3]], [pad_value_at_end_of_vector]]);
    return [res[0][0], res[1][0], res[2][0]];
}

export function unroll_matrix_to_list(matrix) {
    if (!Array.isArray(matrix[0])) {
        return matrix;
    }

    let unrolledArray = [];
    for (let i = 0; i < matrix.length; i++) {
        unrolledArray = unrolledArray.concat(matrix[i]);
    }

    return unrolledArray;
}

