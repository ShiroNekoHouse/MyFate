/*
 * Copyright 2019 The FATE Authors. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.webank.ai.fate.serving.utils;

import com.webank.ai.fate.core.bean.FederatedParty;
import com.webank.ai.fate.core.bean.FederatedRoles;
import com.webank.ai.fate.core.utils.FederatedUtils;
import com.webank.ai.fate.core.utils.ObjectTransform;
import com.webank.ai.fate.serving.core.bean.FederatedInferenceType;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Base64;
import java.util.Map;
import java.util.UUID;

public class InferenceUtils {
    private static final Logger LOGGER = LogManager.getLogger();
    private static final Logger inferenceAuditLogger = LogManager.getLogger("inferenceAudit");
    private static final Logger inferenceLogger = LogManager.getLogger("inference");

    public static String generateCaseid() {
        return UUID.randomUUID().toString().replace("-", "");
    }

    public static String generateSeqno() {
        return UUID.randomUUID().toString().replace("-", "");
    }

    public static void logInference(Enum<FederatedInferenceType> inferenceType, FederatedParty federatedParty, FederatedRoles federatedRoles, String caseid, String seqno, int statusCode, Map<String, Object> featureData) {
        String featureDataBase64String = Base64.getEncoder().encodeToString(ObjectTransform.bean2Json(featureData).getBytes());
        inferenceLogger.info(" {} {} {} {} {} {} {} {}", inferenceType, federatedParty.getRole(), federatedParty.getPartyId(), FederatedUtils.federatedRolesIdentificationString(federatedRoles), caseid, seqno, statusCode, featureDataBase64String);
    }

    public static void logInferenceAudited(Enum<FederatedInferenceType> inferenceType, FederatedParty federatedParty, FederatedRoles federatedRoles, String caseid, String seqno, int statusCode, boolean useCache, boolean billing) {
        String billingDesc, useCacheDesc;
        if (billing) {
            billingDesc = "billing";
        } else {
            billingDesc = "nobilling";
        }
        if (useCache) {
            useCacheDesc = "fromCache";
        } else {
            useCacheDesc = "fromRequest";
        }
        inferenceAuditLogger.info(" {} {} {} {} {} {} {} {} {}", inferenceType, federatedParty.getRole(), federatedParty.getPartyId(), FederatedUtils.federatedRolesIdentificationString(federatedRoles), caseid, seqno, statusCode, useCacheDesc, billingDesc);
    }

    public static Object getClassByName(String classPath) {
        try {
            Class thisClass = Class.forName(classPath);
            return thisClass.getConstructor().newInstance();
        } catch (ClassNotFoundException ex) {
            LOGGER.error("Can not found this class: {}.", classPath);
        } catch (NoSuchMethodException ex) {
            LOGGER.error("Can not get this class({}) constructor.", classPath);
        } catch (Exception ex) {
            LOGGER.error(ex);
            LOGGER.error("Can not create class({}) instance.", classPath);
        }
        return null;
    }
}
