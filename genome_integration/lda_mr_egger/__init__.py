import numpy as np
import scipy.stats as stats
from genome_integration import ivw

class LDAMREgger(ivw.IVWResult):

    def do_lda_mr_egger(self, ld_matrix):

        return self.do_lda_mr_egger_on_estimates(self.outcome_tuples, self.exposure_tuples, ld_matrix)



    def do_lda_mr_egger_on_estimates(self, list_of_outcome_tuples, list_of_exposure_tuples, pearson_ld_matrix, write_out=False):
        """
        This will do LDA mr egger regression as described in Barfield et al. 2018
        Implemented based on their paper, and a reference implementation they provided in personal communication

        <Begin email.>
        Hi Adriaan,
        Below please find the R function to implement the approach. Please let me know if you have any questions.

        -Richard

        X is the vector of  joint eQTL effects
        Y is the vector of joint GWAS effects
        W is the inverse of the covariance of the joint GWAS effects (i.e. var(Y))

        weight.func2<-function(X,Y,W){
          bX<-cbind(1,X)
          bread<-solve(crossprod(bX,W)%*%bX)
          theEsts<-bread%*%crossprod(bX,W%*%Y)
          theresid<-c(Y-theEsts[1]-X*theEsts[2])
          Sig.Est<-c(crossprod(theresid,W%*%theresid))/(length(X)-2)
          finresults<- cbind(theEsts,diag(bread)*Sig.Est)
          TestStat<-theEsts/sqrt(finresults[,2])
          Pvals<-2*pt(abs(TestStat),df = nrow(bX)-2,lower.tail = F)
          return(cbind(finresults,TestStat,Pvals))
        }
        <End email.>


        :param list_of_outcome_tuples:
        :param list_of_exposure_tuples:
        :param pearson_ld_matrix:
        :return:
        """


        if len(list_of_outcome_tuples) < 3:
            raise ValueError("Could not do lda mr egger on estimates, too little estimates supplied")

        marginal_exposure = np.asarray(list_of_exposure_tuples)
        marginal_outcome = np.asarray(list_of_outcome_tuples)

        #flip to make exposure strictly positive.
        to_flip = marginal_exposure[:, 0] < 0

        marginal_exposure[to_flip, 0] = marginal_exposure[to_flip,0] * -1
        marginal_outcome[to_flip, 0] = marginal_outcome[to_flip, 0] * -1
        pearson_ld_matrix[to_flip,:][:,to_flip] = pearson_ld_matrix[to_flip,:][:,to_flip] * -1


        if marginal_exposure.shape[1] < 1:
            raise ValueError("No standard errors supplied to the marginal exposure")

        if marginal_outcome.shape[1] < 1:
            raise ValueError("No standard errors supplied to the marginal outcome")

        sigma = pearson_ld_matrix
        inv_sigma = np.linalg.inv(sigma)

        conditional_outcome = inv_sigma @ marginal_outcome[:,0]

        conditional_exposure = inv_sigma @ marginal_exposure[:,0]

        sigma_g = inv_sigma * np.median(marginal_outcome[:,1])

        b_x = np.concatenate((np.ones((conditional_exposure.shape[0],1)) , conditional_exposure.reshape(conditional_exposure.shape[0], 1)), axis=1)

        bread  = np.linalg.inv(b_x.transpose() @ sigma_g @ b_x)

        estimates = bread @ b_x.transpose() @ sigma_g @ conditional_outcome

        residuals = conditional_outcome - estimates[0] - conditional_exposure * estimates[1]

        significant_estimates = (residuals.transpose() @ (sigma_g @ residuals)) / (conditional_exposure.shape[0] - 2)

        test_stat = estimates / np.sqrt(significant_estimates * np.diag(bread))

        p_val = 2 * stats.t.sf(np.abs(test_stat), df=conditional_exposure.shape[0]-2)

        """
        This only used to compare to the R implementation.
        """
        if write_out:
            with open("check_r_implementation.txt", "w") as f:
                for i in range(marginal_outcome.shape[0]):
                    f.write("{}\t{}\t".format(conditional_exposure[i], conditional_outcome[i]) + "\t".join([str(x) for x in sigma_g[i,:]]) + "\n" )

        return estimates, np.sqrt(np.diag(bread) * significant_estimates), test_stat, p_val

