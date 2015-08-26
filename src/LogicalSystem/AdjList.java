package LogicalSystem;

import java.util.HashMap;
import java.util.HashSet;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 16/08/13
 * Time: 3:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class AdjList {
	private HashMap<Integer, HashSet<Integer>> forwardsAdjList;
	private HashMap<Integer, HashSet<Integer>> reverseAdjList;
	private HashSet<Integer> inputs;
	private HashSet<Integer> outputs;

	public AdjList() {
		this.forwardsAdjList = new HashMap<Integer, HashSet<Integer>>();
		this.reverseAdjList = new HashMap<Integer, HashSet<Integer>>();
		this.inputs = new HashSet<Integer>();
		this.outputs = new HashSet<Integer>();
	}

	public void addEdge(int source, int target) {
		//Check to see if the forwards adj list map contains source node.
		if (forwardsAdjList.containsKey(source)) {
			//If it does, just add another edge.
			forwardsAdjList.get(source).add(target);
		} else {
			//Otherwise, source is a new source node - add it to the inputs set and a new integer
			//set to the adj list.
			HashSet<Integer> iset = new HashSet<Integer>();
			iset.add(target);
			forwardsAdjList.put(source, iset);
		}
		//if the sink list doesn't contain source, it's an input at this point.
		if (!reverseAdjList.containsKey(source) && !reverseAdjList.containsKey(Math.abs(source))) {
			inputs.add(source);
		}
		//Add the reverse edge in.
		if (reverseAdjList.containsKey(target)) {
			reverseAdjList.get(target).add(source);
		} else {
			HashSet<Integer> iset = new HashSet<Integer>();
			iset.add(source);
			reverseAdjList.put(target, iset);
		}
		//Check to see if the forwards adj list map contains the target. If it doesn't, the target is an output.
		if (!forwardsAdjList.containsKey(target) && !forwardsAdjList.containsKey(Math.abs(target))) {
			outputs.add(target);
		}
		//The target node can't be an input, and the source can't be an output.
		inputs.remove(target);
		inputs.remove(Math.abs(target));
		outputs.remove(source);
		outputs.remove(Math.abs(target));
	}

	public HashSet<Integer> getInputs() {
		return inputs;
	}

	public HashSet<Integer> getOutputs() {
		return outputs;
	}

	public boolean isInput(int index) {
		return inputs.contains(index) || inputs.contains(Math.abs(index));
	}

	public boolean isOutput(int index) {
		return outputs.contains(index) || outputs.contains(Math.abs(index));
	}

	public void removeEdge(int source, int target) {
		if (forwardsAdjList.containsKey(source)) {
			forwardsAdjList.get(source).remove(target);
			if (forwardsAdjList.get(source).isEmpty()) {
				forwardsAdjList.remove(source);
			}
		}
		if (reverseAdjList.containsKey(target)) {
			reverseAdjList.get(target).remove(source);
			if (reverseAdjList.get(target).isEmpty()) {
				reverseAdjList.remove(target);
			}
		}
	}

	public HashSet<Integer> getNeighbours(int index) {
		HashSet<Integer> iset = new HashSet<Integer>();
		if (forwardsAdjList.containsKey(index)) {
			iset.addAll(forwardsAdjList.get(index));
		}
		if (forwardsAdjList.containsKey(Math.abs(index))) {
			iset.addAll(forwardsAdjList.get(index));
		}
		if (reverseAdjList.containsKey(index)) {
			iset.addAll(reverseAdjList.get(index));
		}
		if (reverseAdjList.containsKey(Math.abs(index))) {
			iset.addAll(reverseAdjList.get(index));
		}
		return iset;
	}

	public HashSet<Integer> getChildren(int index) {
		HashSet<Integer> returns = new HashSet<Integer>();
		if (forwardsAdjList.containsKey(index)) {
			returns.addAll(forwardsAdjList.get(index));
		}
		if (forwardsAdjList.containsKey(Math.abs(index))) {
			returns.addAll(forwardsAdjList.get(index));
		}
		return returns;
	}

	public HashSet<Integer> getParents(int index) {
		HashSet<Integer> returns = new HashSet<Integer>();
		if (reverseAdjList.containsKey(index)) {
			returns.addAll(reverseAdjList.get(index));
		}
		if (reverseAdjList.containsKey(Math.abs(index))) {
			returns.addAll(reverseAdjList.get(index));
		}
		return returns;
	}

}
