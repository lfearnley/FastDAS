package LogicalSystem;


import org.apache.log4j.Logger;
import org.gnu.glpk.GLPK;
import org.gnu.glpk.GLPKConstants;
import org.gnu.glpk.glp_prob;
import org.gnu.glpk.glp_smcp;

import java.util.concurrent.Callable;

/**
 * A Callable intended for use with single LP instances.
 */
public class LpSolveCallable implements Callable {

	private glp_prob lp;
	private org.apache.log4j.Logger log = Logger.getLogger(LpSolveCallable.class);

	public LpSolveCallable(glp_prob lp) {
		this.lp = lp;
	}

	public double[] call() {
		double[] returnarr = new double[1];
		try {
			glp_smcp parm = new glp_smcp();
			GLPK.glp_init_smcp(parm);
			GLPK.glp_java_set_msg_lvl(GLPKConstants.GLP_MSG_OFF);
			GLPK.glp_write_lp(lp, null, "out.lp");
			int ret = GLPK.glp_simplex(lp, parm);
			int status = GLPK.glp_get_status(lp);
			returnarr = new double[GLPK.glp_get_num_cols(lp) + 1];
			if (ret == 0 && status == GLPK.GLP_OPT) {
				for (int i = 1; i < returnarr.length; i++) {
					returnarr[i] = GLPK.glp_get_col_prim(lp, i);
					//System.out.println(i + GLPK.glp_get_col_name(lp, i) + "\t" + GLPK.glp_get_col_prim(lp, i));
				}
			} else if (ret == 0 && status != GLPK.GLP_OPT) {
				System.exit(0);
				System.out.println("NON-OPTIMAL SOLUTION, RET: " + ret);
				returnarr[0] = -1.0;
			}
		} catch (Exception e) {
			log.error("ERROR : Exception in solver", e);
			System.exit(-1);
		} finally {
			this.lp = null;
			return returnarr;
		}
	}
}
