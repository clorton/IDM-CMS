﻿/***************************************************************************************************

Copyright (c) 2018 Intellectual Ventures Property Holdings, LLC (IVPH) All rights reserved.

EMOD is licensed under the Creative Commons Attribution-Noncommercial-ShareAlike 4.0 License.
To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode

***************************************************************************************************/

using System;
using System.Collections.Generic;
using compartments.emod.interfaces;
using compartments.emod.utils;
using distlib;
using distlib.samplers;

namespace compartments.emod.distributions
{
    public class Normal : INumericOperator, IValue
    {
        private readonly double _mean;
        private readonly double _variance;

        private readonly DistributionSampler _distributionSampler;

        public Normal(double mean, double variance)
        {
            if (variance <= 0.0)
                throw new ArgumentException("Variance must be greater than 0.");

            _mean     = mean;
            _variance = variance;

            _distributionSampler = RandLibSampler.CreateRandLibSampler(RNGFactory.GetRNG());
        }

        public double Value { get { return _distributionSampler.GenerateNormal(_mean, _variance); } }

        public IValue ResolveReferences(IDictionary<string, IValue> map)
        {
            return this;
        }

        public override string ToString()
        {
            return string.Format("(gaussian {0} {1})", _mean, _variance);
        }
    }
}
