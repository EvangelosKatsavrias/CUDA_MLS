function set_polynomialBasis(obj, newPolynomialBasis)
    if (size(obj.m_samplePoints, 2) == newPolynomialBasis.m_dimensions)
        obj.m_polynomial = newPolynomialBasis;
    else error('The number of the basis dimensions does not correspond to the number of the sample points'' coordinates.');
    end
end