/***************************************************************************************************

Copyright (c) 2018 Intellectual Ventures Property Holdings, LLC (IVPH) All rights reserved.

EMOD is licensed under the Creative Commons Attribution-Noncommercial-ShareAlike 4.0 License.
To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode

***************************************************************************************************/

using System;
using System.Collections.Generic;
using System.Linq;
using compartments.emod;
using compartments.solvers.solverbase;
using distlib;
using distlib.samplers;

namespace compartments.solvers
{
    public class BLeaping : SolverBase
    {
        private readonly double[] _currentRates;
        private readonly int[] _firings;
        private readonly double _tau;
        private double _deltaTau;
        private int _step;
        private readonly DistributionSampler _distributionSampler;
        private Dictionary<Species, int> _requires;

        public BLeaping(ModelInfo modelInfo, double duration, int repeats, int samples)
            : base(modelInfo, duration, repeats, samples)
        {
            _currentRates = new double[model.Reactions.Count];
            _firings = new int[model.Reactions.Count];
            _tau = Configuration.CurrentConfiguration.GetParameterWithDefault("b-leaping.Tau", 0.1);
            _step = 1;
            _deltaTau = 0.0;
            _distributionSampler = RandLibSampler.CreateRandLibSampler(rng);
            _requires = new Dictionary<Species, int>();
        }

        protected override void StartRealization()
        {
            base.StartRealization();
            _step = 1;
        }

        protected override double CalculateProposedTau(double tauLimit)
        {
            double desiredTau = _step * _tau;               // Would like to get to current step * step size.
            desiredTau = Math.Min(desiredTau, tauLimit);    // However, don't go farther than the limit
            _deltaTau = desiredTau - CurrentTime;           // Actual step is delta from current time

            UpdateAndSumRates(model.Reactions, _currentRates);

            bool mayProceed = false;
            do
            {
                mayProceed = true;
                /*
                 * for each reaction, calculate a number of firings (Poisson)
                 * if (#firings > 0)
                 *     add #firings to number of reactants needed (by reactant)
                 *     if #reactants needed > #reactants available
                 *         reduce Δt and try again
                */
                foreach (Species s in model.Species) _requires[s] = 0;
                for (int iReaction = 0; iReaction < model.Reactions.Count; ++iReaction)
                {
                    int howMany = _firings[iReaction] = _distributionSampler.GeneratePoisson(_currentRates[iReaction] * _deltaTau);
                    if (howMany > 0)
                    {
                        foreach (Species s in model.Reactions[iReaction].Reactants)
                        {
                            int required = _requires[s] + howMany;
                            if (required > s.Count)
                            {
                                mayProceed = false;
                                break;  // Quit calculating required reactants, we know we have to retry.
                            }
                            else
                            {
                                _requires[s] = required;
                            }
                        }
                        if (!mayProceed)
                        {
                            _deltaTau /= 2.0;
                            break;  // Quit calculating firings, we know we have to retry.
                        }
                    }
                }
            } while (!mayProceed);

            double actualTau = CurrentTime + _deltaTau;

            if (actualTau >= (_step * _tau))
            {
                ++_step;
            }

            return actualTau;
        }

        protected override void ExecuteReactions()
        {
            if (_deltaTau > 0.0)
            {
                for (int jReaction = 0; jReaction < model.Reactions.Count; jReaction++)
                {
                    FireReaction(model.Reactions[jReaction], _firings[jReaction]);
                }
            }
        }

        public override string ToString()
        {
            return "B-Leaping";
        }
    }
}
