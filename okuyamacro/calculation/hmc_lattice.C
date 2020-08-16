double action(void)
{
	int i, V;
	double J;

	int S=0;
	for(i=0;i<V;i++)//loop over all sites
	{
		//Sum over neighbors in positive direction
		J=0.0;
		for(int mu=0;mu<D;mu++) J+=phi
