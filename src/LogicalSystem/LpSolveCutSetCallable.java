package LogicalSystem;

import org.apache.log4j.Logger;

import java.util.HashSet;
import java.util.concurrent.Callable;

import org.gnu.glpk.GLPK;
import org.gnu.glpk.glp_prob;
import org.gnu.glpk.glp_smcp;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 21/08/13
 * Time: 5:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class LpSolveCutSetCallable implements Callable {

	private glp_prob lp;
	private org.apache.log4j.Logger log = Logger.getLogger(LpSolveCallable.class);
	private HashSet<Integer> cutSet = new HashSet<Integer>();

	public LpSolveCutSetCallable(glp_prob lp, HashSet<Integer> cutSet) {
		this.lp = lp;
		this.cutSet = cutSet;
	}

	public HashSet<Integer> getCutSet() {
		return cutSet;
	}

	public Tuple2<double[], HashSet<Integer>> call() {
		double[] returnarr = new double[1];
		Tuple2<double[], HashSet<Integer>> returnTuple = new Tuple2<double[], HashSet<Integer>>(returnarr, this
				.cutSet);
		try {
			glp_smcp parm = new glp_smcp();
			GLPK.glp_init_smcp(parm);
			int ret = GLPK.glp_simplex(lp, parm);
			if (ret == 0) {
				returnarr = new double[GLPK.glp_get_num_cols(lp)+1];
				if (ret == 0) {
					for (int i = 1; i <= returnarr.length; i++) {
						returnarr[i] = GLPK.glp_get_col_prim(lp, i);
					}
				}
			} else {
				returnarr[0] = -1.0;
			}
		} catch (Exception e) {
			log.error("ERROR : Exception in solver", e);
		} finally {
			return returnTuple;
		}
	}
}
