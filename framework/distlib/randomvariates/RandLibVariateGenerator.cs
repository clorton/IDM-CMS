﻿/***************************************************************************************************

Copyright (c) 2018 Intellectual Ventures Property Holdings, LLC (IVPH) All rights reserved.

EMOD is licensed under the Creative Commons Attribution-Noncommercial-ShareAlike 4.0 License.
To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode

***************************************************************************************************/

using distlib.hidden;

namespace distlib.randomvariates
{
    public class RandLibVariateGenerator : RandomVariateGenerator
    {
        private readonly RandLib _randLib;

        // 1/16,777,216 = 1/0x01000000
        private const float Two24Inv = 0.000000059604644775390625f;
        private const double Two52Inv = 0.00000000000000022204460493;

        private RandLibVariateGenerator(int s1 = 1234567890, int s2 = 123456789)
        {
            _randLib = new RandLib(s1, s2);
        }

        public static RandomVariateGenerator CreateRandLibVariateGenerator(uint[] seedData = null)
        {
            int s1 = 1234567890;
            int s2 = 123456789;

            if (seedData != null)
            {
                if (seedData.Length > 0)
                {
                    s1 = (int)seedData[0];
                    if (seedData.Length > 1)
                    {
                        s2 = (int)seedData[1];
                    }
                }
            }

            var generator = new RandLibVariateGenerator(s1, s2);

            return generator;
        }

        public double GenerateUniformOO()
        {
            int rand;

            do
            {
                rand = _randLib.Random();
            } while (rand > 0x3FFFFFBF);    // 0 < rand <= 0x3FFF FFBF

            rand >>= 6;                     // 0 <= rand <= 0x00FF FFFE
            rand += 1;                      // 0 < rand <= 0x00FF FFFF

            // 0 < rand/0x0100 0000 <= 0x00FF FFFF/0x0100 0000 < 1
            // return ((float)rand)/0x01000000;
            return (double)(rand * Two24Inv);
        }

        public double GenerateUniformOC()
        {
            int rand;

            do
            {
                rand = _randLib.Random();
            } while (rand > 0x3FFFFFFF);    // 0 < rand <= 0x3FFF FFFF

            rand >>= 6;                     // 0 <= rand <= 0x00FF FFFF
            rand += 1;                      // 0 < rand <= 0x0100 0000

            // 0 < rand/0x0100 0000 <= 0x0100 0000/0x0100 0000 <= 1
            // return ((float)(rand)) / 0x01000000;
            return (double)(rand * Two24Inv);
        }

        public double GenerateUniformCO()
        {
            int rand;

            do
            {
                rand = _randLib.Random();
            } while (rand > 0x3FFFFFFF);    // 0 < rand <= 0x3FFF FFC0

            rand >>= 6;                     // 0 <= rand <= 0x00FF FFFF

            // 0 <= rand/0x0100 0000 <= 0x00FF FFFF/0x0100 0000 < 1
            // return ((float)(rand))/0x01000000;
            return (double)(rand * Two24Inv);
        }

        public double GenerateUniformCC()
        {
            int rand;

            do
            {
                rand = _randLib.Random();
            } while (rand > 0x4000003F);    // 0 < rand <= 0x4000 0000

            rand >>= 6;                     // 0 <= rand <= 0x0100 0000

            // 0 <= rand/0x0100 0000 <= 0x0100 0000/0x0100 0000 <= 1
            // return ((float)(rand)) / 0x01000000;
            return (double)(rand * Two24Inv);
        }

        public override string ToString()
        {
            return string.Format("RandLib Random Variate Generator ({0})", _randLib);
        }
    }
}
